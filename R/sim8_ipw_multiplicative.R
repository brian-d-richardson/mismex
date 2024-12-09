#' Simulation 8: IPW with multiplicative measurement error
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
sim8.ipw.multiplicative <- function(
    n,
    rel,
    n.supp,
    k,
    B,
    seed) {

  # for troubleshooting -----------------------------------------------------

  #library(MASS); library(devtools); load_all()
  #n = 800; rel = 0.838; B = 80; seed = 1;

  # define parameters -------------------------------------------------------

  gg <- c(0, 1, 1, 1)                       # Y|A,L parameters
  ipw.formula <- "~A1 + A2"                      # MSM formula
  ps.formula <- "~L + I(L^2)"                    # PS model formula
  inv.link <- inv.ident;                         # MSM link function
  d.inv.link <- d.inv.ident;                     # MSM derivative of link
  mc.seed <- 123                                 # MCCS seed value
  var.e <- 1.0336*(rel^(-1) - 1) / (1.0336 + 0.13^2)
  cov.e <- var.e * diag(2)
  coef.a.l <- matrix(                            # coefs in A|L model
    data = c(0, 0,
             0, 0,
             1, -1),
    nrow = 3, byrow = T)
  var.a.l <- c(1, 1)                       # variance of A|L

  # simulate data -----------------------------------------------------------

  set.seed(seed)                                 # seed for reproducibility
  L <- rnorm(n, 0, 0.36)                         # confounder
  A <- mvrnorm(n = n,                            # true exposure
               mu = c(0, 0),
               Sigma = diag(var.a.l)) +
    cbind(1, L, L^2) %*% coef.a.l
  colnames(A) = paste0("A", 1:2)
  Astar <- A * mvrnorm(n = n,                    # mismeasured exposure
                       m = c(1, 1),
                       Sigma = cov.e)
  Y_mean <- cbind(1, A, L) %*% gg                # mean of outcome Y
  Y <- rnorm(n, 0, 1) + Y_mean                   # outcome Y
  colnames(A) <- colnames(Astar) <- c("A1", "A2")
  dat0 <- data.frame(Y, A, L)                    # oracle data
  datstar <- data.frame(Y, Astar, L)             # mismeasured data

  #apply(A, 2, var) / apply(Astar, 2, var)

  # independent replicate exposure measurements ------------------------------

  Asupp <- head(A, n.supp)
  Astarsupp <- lapply(
    X = 1:n.supp,
    FUN = function(ii) {
      do.call(rbind, replicate(k, Asupp[ii,], simplify = F)) +
        mvrnorm(n = k,
                m = c(0, 0),
                Sigma = cov.e)
    })

  ## estimate measurement error covariance
  Astarsup.mean <- lapply(
    X = Astarsupp,
    FUN = colMeans)
  cov.e.hat <- Reduce("+", lapply(
    X = 1:n.supp,
    FUN = function(ii) {
      xx <- t(Astarsupp[[ii]]) - Astarsup.mean[[ii]]
      xx %*% t(xx)
    })) /
    (n.supp * (k - 1))

  # store values for estimation ---------------------------------------------

  len.A <- ncol(A)                               # dimension of A
  mean.a <- colMeans(A)                          # marginal mean of A
  cov.a <- cov(A)*(n-1)/n                        # marginal covariance of A
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
    n, rel, n.supp, k, B, seed,
    res.NI$est[1:3],
    res.CI$est[1:3],
    sqrt(c(
      diag(res.NI$var)[1:3],
      diag(res.CI$var)[1:3],
      diag(res.NI$bc.var)[1:3],
      diag(res.CI$bc.var)[1:3]
    )))

  names(ret) <- c(
    "n", "rel", "n.supp", "k",  "B", "seed",
    apply(tidyr::expand_grid(
      c("ghat", "stde", "bste"),
      c("NI", "CI"),
      1:3), 1, paste, collapse="."))

  # return named numeric vector of length: 24
  return(ret)
}

