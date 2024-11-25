#' run one IPW simulation with a trivariate exposure, binary outcome, and a
#' linear link function, with estimated measurement error covariance
#'
#' @inheritParams sim.gfmla
#' @param k number of measurement error replicates
#' @param n.supp number of observations with replicate exposure measurements
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
sim.ipw.estvar <- function(n,
                           vare,
                           B,
                           k,
                           n.supp,
                           seed) {

  ## for troubleshooting
  #library(MASS); library(devtools); load_all()
  #n = 800; vare = 0.05; B = 20; seed = 1; k = 5; n.supp = 10
  #n = 800; vare = 0.0001; B = 2; seed = 1; k = 5; n.supp = 10

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

  ## independent replicate exposure measurements
  Asupp <- head(A, n.supp)

  Astarsupp <- lapply(
    X = 1:n.supp,
    FUN = function(ii) {
      do.call(rbind, replicate(k, Asupp[ii,], simplify = F)) +
        mvrnorm(n = k,
                m = c(0, 0, 0),
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

  ## MCCS IPW estimator
  res.CI <- fit.ipw.mccs(data = datstar,
                         args = args.ipw,
                         cov.e = cov.e.hat, B = B, mc.seed = mc.seed,
                         mean.a = colMeans(Astar),
                         cov.a = cov(Astar) - cov.e.hat)

  # combine results: estimates and std errors for 4 parameters
  ret <- c(
    n, vare, B, n.supp, k, seed,
    res.CI$est[1:4],
    sqrt(c(diag(res.CI$var)[1:4],
           diag(res.CI$bc.var)[1:4]
  )))

  # return result (numeric vector of length 40)
  names(ret) <- c(
    "n", "vare", "B", "n.supp", "k", "seed",
    apply(tidyr::expand_grid(
      c("ghat", "stde", "bste"),
      c("CI"),
      1:4), 1, paste, collapse="."))

  return(ret)
}
