#' Simulation 5: g-formula with varying reliability
#'
#' @inheritParams fit.gfmla
#'
#' @param n a positive integer, the sample size
#' @param rel a number in [0,1], the measurement error reliability
#' @param B a positive integer, the number of Monte-Carlo replicates for MCCS
#' @param seed a positive integer, the random number seed
#'
#' @return a named numeric vector with the following entries
#' \itemize{
#' \item{n}
#' \item{vare}
#' \item{B}
#' \item{seed}
#' \item{a1, ..., ak}
#' \item{est.xx: estimates of dose response curve using various methods}
#' \item{ste.xx: standard errors using various methods}
#' \item{bcs.xx: bias-corrected standard errors using various methods}
#' }
#'
#' @export
sim5.gfmla.reliability <- function(
    n,
    a,
    rel,
    B,
    seed) {

  # for troubleshooting -----------------------------------------------------

  #library(MASS); library(devtools); library(dplyr); load_all()
  #n = 800; a = -1:2; rel = 0.5; B = 80; seed = 1;

  # define parameters -------------------------------------------------------

  cov.e <- (rel^(-1) - 1) * (1/3)                 # var(epsilon)
  mc.seed <- 123                                  # MC seed
  inv.link <- inv.ident                           # inverse link
  d.inv.link <- d.inv.ident                       # deriv of inv link
  g <- c(0, 0.25, 0.5, -0.5, 1)                   # outcome model parameters
  formula <- "~A + I(A^2) + I(A^3) + L"           # outcome model formula
  args <- list(formula = formula,                 # model fitting arguments
               inv.link = inv.link,
               d.inv.link = d.inv.link)

  # simulate data -----------------------------------------------------------

  set.seed(seed)
  L <- runif(n)                                                  # confounder
  A <- rnorm(n, L, sqrt(0.25))
  EY <- inv.link(model.matrix(as.formula(formula)) %*% g)        # mean of outcome
  Y <- rnorm(n, EY, 0.16)                                          # outcome
  Astar <- A + rnorm(n, 0, sqrt(cov.e))                          # mismeasured A
  dat0 <- data.frame(Y, A, L)                 # oracle data
  datstar <- data.frame(Y, Astar, L)          # mismeasured data
  colnames(dat0) <- colnames(datstar) <- c("Y", "A", "L")

  #var(A) / var(Astar); var(A^2) / var(Astar^2); var(A^3) / var(Astar^3)

  # estimate E{Y(a)} at grid of a -------------------------------------------

  # g-formula
  gfmla.naive <- fit.gfmla(data = datstar, a = a, args = args)

  # corrected g-formula
  gfmla.mccs <- fit.gfmla.mccs(data = datstar, a = a, args = args,
                               cov.e = cov.e, B = B, mc.seed = mc.seed,
                               start = gfmla.naive$est[1:length(g)])

  # extract estimated E{Y(a)} and std error ---------------------------------

  # (i) naive g-formula
  est.NG <- tail(gfmla.naive$est, length(a))
  ste.NG <- sqrt(tail(diag(gfmla.naive$var), length(a)))
  bcs.NG <- sqrt(tail(diag(gfmla.naive$bc.var), length(a)))

  # (ii) corrected g-formula
  est.CG <- tail(gfmla.mccs$est, length(a))
  ste.CG <- sqrt(tail(diag(gfmla.mccs$var), length(a)))
  bcs.CG <- sqrt(tail(diag(gfmla.mccs$bc.var), length(a)))

  # combine results ---------------------------------------------------------

  ret <- c(n, rel, B, seed, a,
           est.NG, est.CG,
           ste.NG, ste.CG,
           bcs.NG, bcs.CG)

  names(ret) <- c(
    "n", "rel", "B", "seed",
    paste0("a", 1:length(a)),
    apply(tidyr::expand_grid(
      c("est", "ste", "bcs"),
      c("NG", "CG"),
      1:length(a)),
      1, paste, collapse="."))

  # return named numeric vector of length: 32
  return(ret)
}
