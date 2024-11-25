#' run one IPW simulation with nonlinear MSM
#'
#' @inheritParams sim.gfmla
#'
#' @return a named numeric vector with the following entries
#' \itemize{
#' \item{n}
#' \item{vare}
#' \item{B}
#' \item{seed}
#' \item{ghat.OL: oracle linear regression estinates}
#' \item{ghat.NL: naive linear regression estinates}
#' \item{ghat.CL: corrected linear regression estinates}
#' \item{ghat.OG: oracle IPW estinates}
#' \item{ghat.NG: naive IPW estinates}
#' \item{ghat.CG: corrected IPW estinates}
#' \item{stde.OL: oracle linear regression standard errors}
#' \item{stde.NL: naive linear regression standard errors}
#' \item{stde.CL: corrected linear regression standard errors}
#' \item{stde.OG: oracle IPW standard errors}
#' \item{stde.NG: naive IPW standard errors}
#' \item{stde.CG: corrected IPW standard errors}
#' }
#'
#' @export
sim.ipw.nonlinearmsm <- function(n,
                                 vare,
                                 B,
                                 seed) {

  ## for troubleshooting
  #library(MASS); library(devtools); load_all()
  #n = 800; vare = 0.5; B = 80; seed = 1;
  #n = 800; vare = 0.0001; B = 2; seed = 1;

  gg <- c(2, 1.5, -1.5, 1, -1, 1, 0.5, -0.5)     # Y|A,L parameters
  ipw.formula <- "~A1 + I(A1^2) + A2 + I(A2^2)"  # MSM formula
  ps.formula <- "~L"                             # PS model formula
  inv.link <- inv.ident;                         # MSM link function
  d.inv.link <- d.inv.ident;                     # MSM derivative of link
  cov.e <- diag(c(vare, vare))                   # measurement error variance
  mc.seed <- 123                                 # MCCS seed value
  coef.a.l <- matrix(
    data = c(0, 0,
             1, -1),                             # coefs in A|L model
    nrow = 2, byrow = T)
  var.a.l <- c(1, 1)                             # variance of A|L

  ## generate data

  set.seed(seed)                                 # seed for reproducibility
  L <- rnorm(n)                                  # confounder
  A <- mvrnorm(n = n,                            # true exposure
               mu = c(0, 0),
               Sigma = diag(var.a.l)) +
    cbind(1, L) %*% coef.a.l
  Astar <- A + mvrnorm(n = n,                    # mismeasured exposure
                       m = c(0, 0),
                       Sigma = cov.e)
  Y_mean <- cbind(1, A, A^2, L, A*L) %*% gg      # mean of outcome Y
  Y <- rnorm(n, 0, 5) + Y_mean                         # outcome Y
  colnames(A) <- colnames(Astar) <- c("A1", "A2")
  dat0 <- data.frame(Y, A, L)                    # oracle data
  datstar <- data.frame(Y, Astar, L)             # mismeasured data

  #plot(L, A[,1]); plot(L, A[,2])
  #plot(Y_mean, Y)
  #plot(A[,1], Y); plot(A[,2], Y);
  #apply(A, 2, var) / apply(Astar, 2, var)

  ## store values for estimation

  len.A <- ncol(A)                               # dimension of A
  mean.a <- colMeans(A)                          # marginal mean of A
  cov.a <- cov(A)                                # marginal covariance of A
  args.ipw <- list(formula = ipw.formula,        # arguments for fitting IPW
                   ps.formula = ps.formula,
                   inv.link = inv.link,
                   d.inv.link = d.inv.link)

  ## estimate E(A | Astar) for regression calibration

  # estimate means and covariances
  E.A <- colMeans(datstar[,c("A1","A2")])            # E(A)                          # E(A)
  E.L <- mean(datstar$L)                             # E(L)
  Sigma.AA <- cov(datstar[,c("A1","A2")]) - cov.e    # Cov(A)
  Sigma.LA <- cov(datstar[,c("A1","A2")], datstar$L) # Cov(A, L)
  Sigma.LL <- var(datstar$L)                         # Cov(L)

  # estimate E(A | Astar, L)
  E.A.AstarL <- E.A + t(
    cbind(Sigma.AA, Sigma.LA) %*%
      solve(rbind(cbind(Sigma.AA + cov.e, Sigma.LA),
                  cbind(t(Sigma.LA), Sigma.LL))) %*%
      t(as.matrix(datstar[,c("A1", "A2", "L")] -
                    do.call(rbind, replicate(n, c(E.A, E.L), simplify = F)))))

  # create data set for regression calibration
  datrc <- data.frame(Y,
                      A1 = E.A.AstarL[,1],
                      A2 = E.A.AstarL[,2],
                      L)

  ## estimate MSM parameters

  # (i) naive IPW estimator
  res.NI <- fit.ipw(data = datstar,
                    args = args.ipw)

  # (ii) oracle IPW estimator
  res.OI <- fit.ipw(data = dat0,
                    args = args.ipw,
                    start = res.NI$est[1:5])

  # (iii) regression calibration IPW estimator
  res.RI <- fit.ipw(data = datrc,
                    args = args.ipw,
                    start = res.NI$est[1:5])

  # (iii) MCCS IPW estimator
  res.CI <- fit.ipw.mccs(data = datstar,
                         args = args.ipw,
                         cov.e = cov.e, B = B, mc.seed = mc.seed,
                         mean.a = colMeans(Astar),
                         cov.a = cov(Astar) - cov.e,
                         start = res.NI$est[1:5])

  # combine results: estimates and std errors for 4 parameters
  ret <- c(
    n, vare, B, seed,
    res.OI$est[1:5], res.NI$est[1:5],
    res.RI$est[1:5], res.CI$est[1:5],
    sqrt(c(
      diag(res.OI$var)[1:5], diag(res.NI$var)[1:5],
      diag(res.RI$var)[1:5], diag(res.CI$var)[1:5],
      diag(res.OI$bc.var)[1:5], diag(res.NI$bc.var)[1:5],
      diag(res.RI$bc.var)[1:5], diag(res.CI$bc.var)[1:5]
    )))

  # return result (numeric vector of length 40)
  names(ret) <- c(
    "n", "vare", "B", "seed",
    apply(tidyr::expand_grid(
      c("ghat", "stde", "bste"),
      c("OI", "NI", "RI", "CI"),
      1:5), 1, paste, collapse="."))

  return(ret)
}
