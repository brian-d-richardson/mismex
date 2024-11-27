#' Simulation 3: IPW with nonlinear PS model
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
sim3.ipw.nonlinearps <- function(
    n,
    vare,
    B,
    seed) {

  # for troubleshooting -----------------------------------------------------

  #library(MASS); library(devtools); load_all()
  #n = 800; vare = 0.05; B = 80; seed = 1;

  # define parameters -------------------------------------------------------

  gg <- c(2, 1.5, -1.5, -2, 1, 1)                # Y|A,L parameters
  ipw.formula <- "~A1 + A2"                      # MSM formula
  ps.formula <- "~L+ I(L^2)"                     # PS model formula
  inv.link <- inv.ident;                         # MSM link function
  d.inv.link <- d.inv.ident;                     # MSM derivative of link
  cov.e <- diag(c(vare, vare))                   # measurement error variance
  mc.seed <- 123                                 # MCCS seed value
  coef.a.l <- matrix(                            # coefs in A|L model
    data = c(0, 0,
             0, 0,
             1, -1),
    nrow = 3, byrow = T)
  var.a.l <- c(0.09, 0.09)                       # variance of A|L

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

  #plot(L, A[,1]); plot(L, A[,2]); apply(A, 2, var) / apply(Astar, 2, var)

  # store values for estimation ---------------------------------------------

  len.A <- ncol(A)                               # dimension of A
  mean.a <- colMeans(A)                          # marginal mean of A
  cov.a <- cov(A)*(n-1)/n                        # marginal covariance of A
  args.ipw <- list(formula = ipw.formula,        # arguments for fitting IPW
                   ps.formula = ps.formula,
                   inv.link = inv.link,
                   d.inv.link = d.inv.link)


  # estimate E(A | Astar) for regression calibration ------------------------

  # estimate means and covariances
  E.A <- colMeans(datstar[,c("A1","A2")])            # E(A)
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


  # estimate MSM parameters -------------------------------------------------

  # (i) naive IPW estimator
  res.NI <- fit.ipw(data = datstar,
                    args = args.ipw)

  # (ii) oracle IPW estimator
  res.OI <- fit.ipw(data = dat0,
                    args = args.ipw,
                    start = res.NI$est[1:3])

  # (iii) regression calibration IPW estimator
  res.RI <- fit.ipw(data = datrc,
                    args = args.ipw,
                    start = res.NI$est[1:3])

  # (iv) MCCS IPW estimator
  res.CI <- fit.ipw.mccs(data = datstar,
                         args = args.ipw,
                         cov.e = cov.e, B = B, mc.seed = mc.seed,
                         mean.a = colMeans(Astar),
                         cov.a = cov(Astar) - cov.e,
                         start = res.NI$est[1:3])

  # (v) SIMEX IPW estimator
  xi <- seq(0, 2, length = 10) # sequence of xi values
  K <- 20 # number of reps per xi value                # number of reps per xi value

  simex.in <- rbind(
    c(xi = 0, k = 1),
    expand.grid(
      xi = tail(xi, -1),
      k = 1:K))
  simex.ests <- t(vapply(
    X = 1:nrow(simex.in),
    FUN.VALUE = numeric(5),
    FUN = function(jj) {
      dat_jk <- datstar
      dat_jk[,c("A1", "A2")] <- dat_jk[,c("A1", "A2")] +
        mvrnorm(n = n, m = c(0, 0), Sigma = cov.e * simex.in$xi[jj])
      res_jk <- fit.ipw(data = dat_jk,
                        args = args.ipw,
                        start = res.NI$est[1:3],
                        return.var = F,
                        return.bcvar = F)
      c(xi = simex.in$xi[jj],
        k = simex.in$k[jj],
        res_jk$est[1:3])
    }
  )) %>%
    as.data.frame() %>%
    pivot_longer(cols = !c(xi, k))

  extrap.model <- lm(                    # extrapolation model
    data = simex.ests,
    formula = value ~ name + xi + I(xi^2))

  res.SI <- predict(                     # extrapolation to xi = -1
    extrap.model,
    newdata = expand.grid(xi = -1, name = unique(simex.ests$name))
  )

  # combine results ---------------------------------------------------------

  ret <- c(
    n, vare, B, seed,
    res.OI$est[1:3], res.NI$est[1:3],
    res.RI$est[1:3], res.SI,
    res.CI$est[1:3],
    sqrt(c(
      diag(res.OI$var)[1:3], diag(res.NI$var)[1:3],
      diag(res.RI$var)[1:3], rep(NA, 3),
      diag(res.CI$var)[1:3],
      diag(res.OI$bc.var)[1:3], diag(res.NI$bc.var)[1:3],
      diag(res.RI$bc.var)[1:3], rep(NA, 3),
      diag(res.CI$bc.var)[1:3]
    )))

  names(ret) <- c(
    "n", "vare", "B", "seed",
    apply(tidyr::expand_grid(
      c("ghat", "stde", "bste"),
      c("OI", "NI", "RI", "SI", "CI"),
      1:3), 1, paste, collapse="."))

  # return named numeric vector of length: 49
  return(ret)
}

