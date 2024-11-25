#' Simulation 6: IPW with varying reliability
#'
#' @inheritParams sim1.gfmla.nonlinear.coarse
#'
#' @return a named numeric vector with the following entries
#' \itemize{
#' \item{n}
#' \item{vare}
#' \item{B}
#' \item{seed}
#' \item{ghat.xx: MSM parameter estimates for various methods}
#' \item{stde.xx: standard errors for various methods}
#' \item{bste.xx: bias-corrected standard errors for various methods}
#' }
#'
#' @export
sim6.ipw.reliability <- function(
    n,
    rel,
    B,
    seed) {

  # for troubleshooting -----------------------------------------------------

  #library(MASS); library(devtools); load_all()
  #n = 800; rel = 0.5; B = 80; seed = 1;

  # define parameters -------------------------------------------------------

  gg <- c(2, 1.5, -1.5, -2, 1, 1)                # Y|A,L parameters
  glm.formula <- "~A1*L + A2*L"                  # Y|A,L model formula
  ipw.formula <- "~A1 + A2"                      # MSM formula
  ps.formula <- "~L+ I(L^2)"                     # PS model formula
  inv.link <- inv.ident;                         # MSM link function
  d.inv.link <- d.inv.ident;                     # MSM derivative of link
  vare <- (rel^(-1) - 1) * 0.284
  cov.e <- diag(c(vare, vare))                   # measurement error variance
  mc.seed <- 123                                 # MCCS seed value
  coef.a.l <- matrix(                            # coefs in A|L model
    data = c(0, 0,
             0, 0,
             1, -1),
    nrow = 3, byrow = T)
  var.a.l <- c(0.25, 0.25)                       # variance of A|L

  # simulate data -----------------------------------------------------------

  set.seed(seed)                                 # seed for reproducibility
  L <- rnorm(n, 0, 0.36)                         # confounder
  A <- mvrnorm(n = n,                            # true exposure
               mu = c(0, 0),
               Sigma = diag(var.a.l)) +
    cbind(1, L, L^2) %*% coef.a.l
  colnames(A) = paste0("A", 1:2)
  Astar <- A + mvrnorm(n = n,                    # mismeasured exposure
                       m = c(0, 0),
                       Sigma = cov.e)
  Y_mean <- cbind(1, A, L, A*L) %*% gg           # mean of outcome Y
  Y <- rnorm(n, 0, 1) + Y_mean                  # outcome Y
  colnames(A) <- colnames(Astar) <- c("A1", "A2")
  dat0 <- data.frame(Y, A, L)                    # oracle data
  datstar <- data.frame(Y, Astar, L)             # mismeasured data

  # store values for estimation ---------------------------------------------

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

  # estimate MSM parameters -------------------------------------------------

  # (i) naive IPW estimator
  res.NI <- fit.ipw(data = datstar,
                    args = args.ipw)

  # (ii) MCCS IPW estimator
  res.CI <- fit.ipw.mccs(data = datstar,
                         args = args.ipw,
                         cov.e = cov.e, B = B, mc.seed = mc.seed,
                         mean.a = colMeans(Astar),
                         cov.a = cov(Astar) - cov.e,
                         start = res.NI$est[1:3])

  # combine results ---------------------------------------------------------

  ret <- c(
    n, rel, B, seed,
    res.NI$est[1:3],
    res.CI$est[1:3],
    sqrt(c(
      diag(res.NI$var)[1:3],
      diag(res.CI$var)[1:3],
      diag(res.NI$bc.var)[1:3],
      diag(res.CI$bc.var)[1:3]
    )))

  names(ret) <- c(
    "n", "rel", "B", "seed",
    apply(tidyr::expand_grid(
      c("ghat", "stde", "bste"),
      c("NI", "CI"),
      1:3), 1, paste, collapse="."))

  # return named numeric vector of length: 22
  return(ret)
}

