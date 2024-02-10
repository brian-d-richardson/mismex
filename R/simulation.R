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
sim1 <- function(n = 8000,
                 vare = 0.05,
                 B = 50,
                 seed = 1) {

  ## for troubleshooting
  #library(MASS); library(devtools); load_all()
  #n = 8000; vare = 0.05; B = 50; seed = 1;

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
  #range(Y_prob); hist(Y_prob)
  Y_prob[Y_prob < 0] <- 0                        # correct Y_prob in rare cases
  Y_prob[Y_prob > 1] <- 1
  Y <- rbinom(n, 1, Y_prob)                      # binary outcome

  len.a <- ncol(A)                               # dimension of A
  len.L <- 8                                     # length of GLM estimate
  len.I <- 13                                    # length of IPW estimate

  ## estimate MSM parameters

  # (i) oracle logistic regression
  ghat.OL <- tryCatch(
    expr = fit.glm(Y = Y, A = A, L = L,
                   inv.link = inv.link, d.inv.link = d.inv.link),
    warning = function(w) rep(NA, len.L),
    error = function(e) rep(NA, len.L))

  evar.OL <- tryCatch(
    expr = get.sand.est(
      ghat = ghat.OL,
      len.a = len.a, n = n,
      get.psi = function(x) get.psi.glm(
        Y = Y, A = A, L = L, g = x,
        inv.link = inv.link, d.inv.link = d.inv.link, return.sums = F)),
    warning = function(w) matrix(NA, len.L, len.L),
    error = function(e) matrix(NA, len.L, len.L))

  # (ii) naive logistic regression
  ghat.NL <- tryCatch(
    expr = fit.glm(Y = Y, A = Astar, L = L,
                   inv.link = inv.link, d.inv.link = d.inv.link),
    warning = function(w) rep(NA, len.L),
    error = function(e) rep(NA, len.L))

  evar.NL <- tryCatch(
    expr = get.sand.est(
      ghat = ghat.NL,
      len.a = len.a, n = n,
      get.psi = function(x) get.psi.glm(
        Y = Y, A = Astar, L = L, g = x,
        inv.link = inv.link, d.inv.link = d.inv.link, return.sums = F)),
    warning = function(w) matrix(NA, len.L, len.L),
    error = function(e) matrix(NA, len.L, len.L))

  # (iii) MCCS logistic regression
  ghat.CL <- tryCatch(
    expr = fit.glm.mccs(Y = Y, Astar = Astar, L = L, var.e = var.e,
                        inv.link = inv.link, d.inv.link = d.inv.link,
                        B = B, seed = 123),
    warning = function(w) rep(NA, len.L),
    error = function(e) rep(NA, len.L))

  evar.CL <- tryCatch(
    expr = get.sand.est(
      ghat = ghat.CL,
      len.a = len.a, n = n,
      get.psi = function(x) get.psi.glm.mccs(
        Y = Y, Astar = Astar, L = L, g = x,
        var.e = var.e, B = B, seed = 123,
        inv.link = inv.link, d.inv.link = d.inv.link, return.sums = F)),
    warning = function(w) matrix(NA, len.L, len.L),
    error = function(e) matrix(NA, len.L, len.L))

  # (iv) oracle IPW estimator
  ghat.OI <- tryCatch(
    expr = fit.ipw(Y = Y, A = A, L = L,
                   inv.link = inv.link, d.inv.link = d.inv.link),
    warning = function(w) rep(NA, len.I),
    error = function(e) rep(NA, len.I))

  get.psi.ipw(
    Y = Y, A = A, L = L, g = x[1:(len.a + 1)],
    inv.link = inv.link, d.inv.link = d.inv.link,
    coef.a.l = coef.a.l, var.a.l = var.a.l,
    mean.a = mean.a, cov.a = cov.a)

  evar.OI <- tryCatch(
    expr = get.sand.est(
      ghat = ghat.OI,
      len.a = len.a, n = n,
      get.psi = function(x) {
        coef.a.l <- matrix(x[len.a + 1 + 1:(2 * len.a)], nrow = len.a, byrow = F)
        var.a.l <- exp(x[3 * len.a + 1 + 1:len.a])
        cbind(
          get.psi.ipw(
            Y = Y, A = A, L = L, g = x[1:(len.a + 1)],
            inv.link = inv.link, d.inv.link = d.inv.link,
            coef.a.l = coef.a.l, var.a.l = var.a.l,
            mean.a = mean(A), cov.a = cov(A),
            return.sums = F),
          get.psi.ps(
            A = A, L = L,
            coef.a.l = coef.a.l, var.a.l = var.a.l,
            return.sums = F))
      }),
    warning = function(w) matrix(NA, len.I, len.I),
    error = function(e) matrix(NA, len.I, len.I))
  #round(sqrt(diag(evar.OI)[1:4]), 4)

  # (iv) naive IPW estimator
  ghat.NI <- tryCatch(
    expr = fit.ipw(Y = Y, A = Astar, L = L,
                   inv.link = inv.link, d.inv.link = d.inv.link),
    warning = function(w) rep(NA, len.I),
    error = function(e) rep(NA, len.I))

  evar.NI <- tryCatch(
    expr = get.sand.est(
      ghat = ghat.NI,
      len.a = len.a, n = n,
      get.psi = function(x) {
        cbind(
          get.psi.ipw(
            Y = Y, A = Astar, L = L, g = x[1:((len.a) + 1)],
            inv.link = inv.link, d.inv.link = d.inv.link,
            coef.a.l = matrix(x[len.a + 1 + 1:(2 * len.a)], nrow = len.a, byrow = F),
            var.a.l = exp(x[3 * len.a + 1 + 1:len.a]),
            mean.a = mean(A), cov.a = cov(A),
            return.sums = F),
          get.psi.ps(
            A = Astar, L = L,
            coef.a.l = matrix(x[len.a + 1 + 1:(2 * len.a)], nrow = len.a, byrow = F),
            var.a.l = exp(x[3 * len.a + 1 + 1:len.a]),
            return.sums = F))
      }),
    warning = function(w) matrix(NA, len.I, len.I),
    error = function(e) matrix(NA, len.I, len.I))


  # (vi) MCCS IPW estimator
  ghat.CI <- tryCatch(
    expr = fit.ipw.mccs(Y = Y, Astar = Astar, L = L,
                        var.e = var.e, B = B, seed = 123,
                        inv.link = inv.link, d.inv.link = d.inv.link),
    warning = function(w) rep(NA, len.I),
    error = function(e) rep(NA, len.I))

  evar.CI <- tryCatch(
    expr = get.sand.est(
      ghat = ghat.CI,
      len.a = len.a, n = n,
      get.psi = function(x) {
        coef.a.l <- matrix(x[len.a + 1 + 1:(2 * len.a)], nrow = len.a, byrow = F)
        var.a.l <- exp(x[3 * len.a + 1 + 1:len.a])
        cbind(
          get.psi.ipw.mccs(
            Y = Y, Astar = Astar, L = L, g = x[1:((len.a) + 1)],
            var.e = var.e, B = B, seed = 123,
            inv.link = inv.link, d.inv.link = d.inv.link,
            coef.a.l = coef.a.l,
            var.a.l = var.a.l,
            mean.a = mean(Astar), cov.a = cov(Astar) - diag(var.e),
            return.sums = F),
          get.psi.ps(
            A = Astar, L = L,
            coef.a.l = coef.a.l,
            var.a.l = var.a.l + var.e,
            return.sums = F))
      }),
    warning = function(w) matrix(NA, len.I, len.I),
    error = function(e) matrix(NA, len.I, len.I))

  ret <- c(n, vare, B, seed,
           ghat.OL[1:4], ghat.NL[1:4], ghat.CL[1:4],
           ghat.OI[1:4], ghat.NI[1:4], ghat.CI[1:4],
           sqrt(c(
             diag(evar.OL)[1:4], diag(evar.NL)[1:4], diag(evar.CL)[1:4],
             diag(evar.OI)[1:4], diag(evar.NI)[1:4], diag(evar.CI)[1:4]
           )))

  names(ret) <- c(
    "n", "vare", "B", "seed",
    apply(tidyr::expand_grid(
      c("ghat", "stde"),
      c("OL", "NL", "CL", "OI", "NI", "CI"),
      1:4), 1, paste, collapse="."))

  return(ret)
}

#library(tictoc); tic("one sim"); sim.res <- sim2(); toc()
#round(sim.res, 2)




