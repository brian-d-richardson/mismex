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

  ## estimate MSM parameters

  # (i) naive logistic regression
  res.NL <- fit.glm(data = datstar,
                    args = args.glm)

  # (ii) oracle logistic regression
  res.OL <- fit.glm(data = dat0,
                    args = args.glm,
                    start = res.NL$est)

  # (iii) MCCS logistic regression
  res.CL <- fit.glm.mccs(dat = dat0,
                         args = args.glm,
                         cov.e = cov.e, B = B, mc.seed = mc.seed,
                         start = res.NL$est)

  # (iv) naive IPW estimator
  res.NI <- fit.ipw(data = datstar,
                    args = args.ipw,
                    mean.a = mean.a,
                    cov.a = cov.a,
                    start = res.NL$est[1:4])

  # (v) oracle IPW estimator
  res.OI <- fit.ipw(data = dat0,
                    args = args.ipw,
                    mean.a = mean.a,
                    cov.a = cov.a,
                    start = res.NI$est[1:4])

  # (vi) MCCS IPW estimator
  res.CI <- fit.ipw.mccs(data = datstar,
                         args = args.ipw,
                         cov.e = cov.e, B = B, mc.seed = mc.seed,
                         mean.a = colMeans(Astar),
                         cov.a = cov(Astar) - cov.e,
                         start = res.NI$est[1:4])

  # combine results: estimates and std errors for 4 parameters
  ret <- c(n, vare, B, seed,
           res.OL$est[1:4], res.NL$est[1:4], res.CL$est[1:4],
           res.OI$est[1:4], res.NI$est[1:4], res.CI$est[1:4],
           sqrt(c(
             diag(res.OL$var)[1:4], diag(res.NL$var)[1:4], diag(res.CL$var)[1:4],
             diag(res.OI$var)[1:4], diag(res.NI$var)[1:4], diag(res.CI$var)[1:4]
           )))

  names(ret) <- c(
    "n", "vare", "B", "seed",
    apply(tidyr::expand_grid(
      c("ghat", "stde"),
      c("OL", "NL", "CL", "OI", "NI", "CI"),
      1:4), 1, paste, collapse="."))

  round(ret, 2)

  return(ret)
}
