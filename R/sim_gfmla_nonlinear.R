#' run one G-formula nonlinear model simulation
#'
#' @inheritParams sim.dr
#' @inheritParams fit.gfmla
#'
#' @param vare a non-negative number, the measurement error variance
#'
#' @return a named numeric vector with the following entries
#' \itemize{
#' \item{n}
#' \item{vare}
#' \item{B}
#' \item{seed}
#' \item{a1, ..., ak}
#' \item{est.OL: oracle logistic regression estinates}
#' \item{est.NL: naive logistic regression estinates}
#' \item{est.CL: corrected logistic regression estinates}
#' \item{est.OG: oracle g-formula estinates}
#' \item{est.NG: naive g-formula estinates}
#' \item{est.CG: corrected g-formula estinates}
#' \item{ste.OL: oracle linear regression standard errors}
#' \item{ste.NL: naive linear regression standard errors}
#' \item{ste.CL: corrected linear regression standard errors}
#' \item{ste.OG: oracle g-formula standard errors}
#' \item{ste.NG: naive g-formula standard errors}
#' \item{ste.CG: corrected g-formula standard errors}
#' }
#'
#' @export
sim.gfmla.nonlinear <- function(n,
                                a,
                                vare,
                                B,
                                seed) {

  ## for troubleshooting
  #library(MASS); library(devtools); load_all()
  #n = 800; a = 3; vare = 0.25; B = 80; seed = 1;
  #n = 800; a = 0:4; vare = 0.0025; B = 2; seed = 1;

  cov.e <- vare                                   # var(epsilon)
  mc.seed <- 123                                  # MC seed
  inv.link <- inv.ident                           # inverse link
  d.inv.link <- d.inv.ident                       # deriv of inv link
  g <- c(0, 0.25, 0.5, -0.5, 1)                   # outcome model parameters
  formula <- "~A + I(A^2) + I(A^3) + L"           # outcome model formula
  args <- list(formula = formula,                 # model fitting arguments
               inv.link = inv.link,
               d.inv.link = d.inv.link)

  # simulate data
  set.seed(seed)
  L <- runif(n)                                                  # confounder
  A <- rnorm(n, L, sqrt(0.25))
  EY <- inv.link(model.matrix(as.formula(formula)) %*% g)        # mean of outcome
  Y <- rnorm(n, EY, 0.16)                                          # outcome
  Astar <- A + rnorm(n, 0, sqrt(cov.e))                          # mismeasured A
  dat0 <- data.frame(Y, A, L)                 # oracle data
  datstar <- data.frame(Y, Astar, L)          # mismeasured data
  colnames(dat0) <- colnames(datstar) <- c("Y", "A", "L")

  # estimate E(A | Astar) for regression calibration ------------------------

  # estimate means and covariances
  E.A <- mean(datstar$A)                  # E(A)
  E.L <- mean(datstar$L)                  # E(L)
  Sigma.AA <- var(datstar$A) - cov.e      # Cov(A)
  Sigma.LA <- cov(datstar$A, datstar$L)   # Cov(A, L)
  Sigma.LL <- var(datstar$L)              # Cov(L)

  # estimate E(A | Astar, L)
  E.A.AstarL <- E.A + c(
    c(Sigma.AA, Sigma.LA) %*%
      solve(rbind(cbind(Sigma.AA + cov.e, Sigma.LA),
                  cbind(t(Sigma.LA), Sigma.LL))) %*%
      t(as.matrix(datstar[,c("A", "L")] -
                    do.call(rbind, replicate(n, c(E.A, E.L), simplify = F)))))

  # create data set for regression calibration
  datrc <- data.frame(Y, A = E.A.AstarL, L)

  # estimate E{Y(a)} at grid of a -------------------------------------------

  # g-formula
  gfmla.naive <- fit.gfmla(data = datstar, a = a, args = args)

  # oracle g-formula
  gfmla.oracle <- fit.gfmla(data = dat0, a = a, args = args,
                            start = gfmla.naive$est[1:length(g)])

  # regression calibration g-formula
  gfmla.rc <- fit.gfmla(data = datrc, a = a, args = args,
                        start = gfmla.naive$est[1:length(g)])

  # corrected g-formula
  gfmla.mccs <- fit.gfmla.mccs(data = datstar, a = a, args = args,
                               cov.e = cov.e, B = B, mc.seed = mc.seed,
                               start = gfmla.naive$est[1:length(g)])

  ## extract estimated E{Y(a)} and std error

  # (i) oracle g-formula
  est.OG <- tail(gfmla.oracle$est, length(a))
  ste.OG <- sqrt(tail(diag(gfmla.oracle$var), length(a)))
  bcs.OG <- sqrt(tail(diag(gfmla.oracle$bc.var), length(a)))

  # (ii) naive g-formula
  est.NG <- tail(gfmla.naive$est, length(a))
  ste.NG <- sqrt(tail(diag(gfmla.naive$var), length(a)))
  bcs.NG <- sqrt(tail(diag(gfmla.naive$bc.var), length(a)))

  # (iii) regression calibration g-formula
  est.RG <- tail(gfmla.rc$est, length(a))
  ste.RG <- sqrt(tail(diag(gfmla.rc$var), length(a)))
  bcs.RG <- sqrt(tail(diag(gfmla.rc$bc.var), length(a)))

  # (iv) corrected g-formula
  est.CG <- tail(gfmla.mccs$est, length(a))
  ste.CG <- sqrt(tail(diag(gfmla.mccs$var), length(a)))
  bcs.CG <- sqrt(tail(diag(gfmla.mccs$bc.var), length(a)))

  # combine results: estimates and std errors for 4 parameters
  ret <- c(n, vare, B, seed, a,
           est.OG, est.NG, est.RG, est.CG,
           ste.OG, ste.NG, ste.RG, ste.CG,
           bcs.OG, bcs.NG, bcs.RG, bcs.CG)

  names(ret) <- c(
    "n", "vare", "B", "seed",
    paste0("a", 1:length(a)),
    apply(tidyr::expand_grid(
      c("est", "ste", "bcs"),
      c("OG", "NG", "RG", "CG"),
      1:length(a)),
      1, paste, collapse="."))

  return(ret)
}
