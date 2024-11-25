#' run one IPW simulation with a trivariate exposure, binary outcome, and a linear link function
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

  ## estimate E(A | Astar) for regression calibration

  # estimate means and covariances
  E.A <- colMeans(datstar[,c("A1","A2","A3")])            # E(A)                          # E(A)
  E.L <- mean(datstar$L)                                  # E(L)
  Sigma.AA <- cov(datstar[,c("A1","A2","A3")]) - cov.e    # Cov(A)
  Sigma.LA <- cov(datstar[,c("A1","A2","A3")], datstar$L) # Cov(A, L)
  Sigma.LL <- var(datstar$L)                              # Cov(L)

  # estimate E(A | Astar, L)
  E.A.AstarL <- E.A + t(
    cbind(Sigma.AA, Sigma.LA) %*%
      solve(rbind(cbind(Sigma.AA + cov.e, Sigma.LA),
                  cbind(t(Sigma.LA), Sigma.LL))) %*%
      t(as.matrix(datstar[,c("A1", "A2", "A3", "L")] -
                  do.call(rbind, replicate(n, c(E.A, E.L), simplify = F)))))

  # create data set for regression calibration
  datrc <- data.frame(Y,
                      A1 = E.A.AstarL[,1],
                      A2 = E.A.AstarL[,2],
                      A3 = E.A.AstarL[,3],
                      L)

  ## estimate MSM parameters

  # (i) naive IPW estimator
  res.NI <- fit.ipw(data = datstar,
                    args = args.ipw)

  # (ii) oracle IPW estimator
  res.OI <- fit.ipw(data = dat0,
                    args = args.ipw,
                    start = res.NI$est[1:4])

  # (iii) regression calibration IPW estimator
  res.RI <- fit.ipw(data = datrc,
                    args = args.ipw,
                    start = res.NI$est[1:4])

  # (iv) MCCS IPW estimator
  res.CI <- fit.ipw.mccs(data = datstar,
                         args = args.ipw,
                         cov.e = cov.e, B = B, mc.seed = mc.seed,
                         mean.a = colMeans(Astar),
                         cov.a = cov(Astar) - cov.e,
                         start = res.NI$est[1:4])

  # (v) SIMEX IPW estimator
  xi <- seq(0, 2, length = 5) # sequence of xi values
  K <- 50                     # number of reps per xi value

  simex.in <- expand.grid(    # esimates of MSM params with simulated error
    xi = xi,
    k = 1:K)
  simex.ests <- t(vapply(
    X = 1:nrow(simex.in),
    FUN.VALUE = numeric(6),
    FUN = function(jj) {
      dat_jk <- datstar
      dat_jk[,c("A1", "A2", "A3")] <- dat_jk[,c("A1", "A2", "A3")] +
        mvrnorm(n = n, m = c(0, 0, 0), Sigma = cov.e * simex.in$xi[jj])
      res_jk <- fit.ipw(data = dat_jk,
                        args = args.ipw,
                        start = res.NI$est[1:4],
                        return.var = F,
                        return.bcvar = F)
      c(xi = simex.in$xi[jj],
        k = simex.in$k[jj],
        res_jk$est[1:4])
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

  # combine results: estimates and std errors for 4 parameters
  ret <- c(
    n, vare, B, seed,
    res.OI$est[1:4], res.NI$est[1:4],
    res.RI$est[1:4], res.SI, res.CI$est[1:4],
    sqrt(c(
    diag(res.OI$var)[1:4], diag(res.NI$var)[1:4],
    diag(res.RI$var)[1:4], rep(NA, 4), diag(res.CI$var)[1:4],
    diag(res.OI$bc.var)[1:4], diag(res.NI$bc.var)[1:4],
    diag(res.RI$bc.var)[1:4], rep(NA, 4), diag(res.CI$bc.var)[1:4]
  )))

  # return result (numeric vector of length 40)
  names(ret) <- c(
    "n", "vare", "B", "seed",
    apply(tidyr::expand_grid(
      c("ghat", "stde", "bste"),
      c("OI", "NI", "RI", "SI", "CI"),
      1:4), 1, paste, collapse="."))

  return(ret)
}
