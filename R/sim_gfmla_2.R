#' run one G-formula 2 simulation
#'
#' @inheritParams sim.dr
#' @inheritParams fit.gfmla
#'
#' @return a named numeric vector with the following entries
#' \itemize{
#' \item{n}
#' \item{vare}
#' \item{B}
#' \item{seed}
#' \item{a1, ..., ak}
#' \item{est.OG: oracle g-formula estinates}
#' \item{est.NG: naive g-formula estinates}
#' \item{est.CG: corrected g-formula estinates}
#' }
#'
#' @export
sim.gfmla.2 <- function(n,
                      a,
                      vare,
                      B,
                      seed) {

  ## for troubleshooting
  #library(MASS); library(devtools); load_all()
  #n = 8000; a = seq(1, 4, length = 20); vare = 0.25; B = 80; seed = 1;
  #n = 8000; a = seq(1, 4, length = 20); vare = 0.0025; B = 2; seed = 1;

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

  # estimate E(A | Astar) for regression calibration ------------------------

  # estimate means and covariances
  E.A <- mean(datstar$A)                               # E(A)
  E.L <- colMeans(datstar[,c("L1","L2")])              # E(L)
  Sigma.AA <- var(datstar$A) - cov.e                   # Cov(A)
  Sigma.LA <- cov(datstar$A, datstar[,c("L1","L2")])   # Cov(A, L)
  Sigma.LL <- cov(datstar[,c("L1","L2")])              # Cov(L)

  # estimate E(A | Astar, L)
  E.A.AstarL <- E.A + c(
    c(Sigma.AA, Sigma.LA) %*%
      solve(rbind(cbind(Sigma.AA + cov.e, Sigma.LA),
                  cbind(t(Sigma.LA), Sigma.LL))) %*%
      t(as.matrix(datstar[,c("A", "L1", "L2")] -
                    do.call(rbind, replicate(n, c(E.A, E.L), simplify = F)))))

  # create data set for regression calibration
  datrc <- data.frame(Y, A = E.A.AstarL, L1, L2)

  # estimate E{Y(a)} at grid of a -------------------------------------------

  # naive g-formula
  gfmla.naive <- fit.gfmla(data = datstar, a = a, args = args,
                           return.var = F, return.bcvar = F)

  # oracle g-formula
  gfmla.oracle <- fit.gfmla(data = dat0, a = a, args = args,
                            start = gfmla.naive$est[1:length(g)],
                            return.var = F, return.bcvar = F)

  # regression calibration g-formula
  gfmla.rc <- fit.gfmla(data = datrc, a = a, args = args,
                        start = gfmla.naive$est[1:length(g)],
                        return.var = F, return.bcvar = F)

  # corrected g-formula
  gfmla.mccs <- fit.gfmla.mccs(data = datstar, a = a, args = args,
                               cov.e = cov.e, B = B, mc.seed = mc.seed,
                               start = gfmla.naive$est[1:length(g)],
                               return.var = F, return.bcvar = F)

  ## extract estimated E{Y(a)}

  # (i) oracle g-formula
  est.OG <- tail(gfmla.oracle$est, length(a))

  # (ii) naive g-formula
  est.NG <- tail(gfmla.naive$est, length(a))

  # (iii) regression calibration g-formula
  est.RG <- tail(gfmla.rc$est, length(a))

  # (iv) corrected g-formula
  est.CG <- tail(gfmla.mccs$est, length(a))

  # combine results: estimates and std errors for 4 parameters
  ret <- c(n, vare, B, seed, a,
           est.OG, est.NG, est.RG, est.CG)

  names(ret) <- c(
    "n", "vare", "B", "seed",
    paste0("a", 1:length(a)),
    apply(tidyr::expand_grid(
      c("OG", "NG", "RG", "CG"),
      1:length(a)),
      1, paste, collapse="."))

  return(ret)
}
