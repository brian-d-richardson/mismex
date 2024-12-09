#' Simulation 10: DR with case cohort sampling
#'
#' @inheritParams sim1.gfmla.nonlinear.coarse
#'
#' @param pi.cc a number in (0,1], the case-cohort sampling proportion
#'
#' @return a data frame with the following columns
#' \itemize{
#' \item{n}
#' \item{B}
#' \item{vare}
#' \item{pi.cc}
#' \item{seed}
#' \item{ps: an indicator for whether the propensity score model is correctly
#' (0) or incorrectly (1) specified}
#' \item{out: an indicator for whether the outcome model is correctly
#' (0) or incorrectly (1) specified}
#' \item{type: the type of estimator (naive, oracle, reg. cal., SIMEX, or MCCS)}
#' \item{est: estimated MSM parameter}
#' \item{se: standard error}
#' \item{bse: bias-corrected standard error}
#' }
#'
#' @export
sim10.dr.cc <- function(
    n,
    vare,
    B,
    pi.cc,
    seed) {

  # for troubleshooting -----------------------------------------------------

  #library(MASS); library(devtools); load_all();
  #n = 2000; vare = 0.02; B = 80; pi.cc = 0.5; seed = 1;

  # define parameters -------------------------------------------------------

  mc.seed <- 123                                # MCCS seed
  inv.link <- inv.ident                         # inverse link
  d.inv.link <- d.inv.ident                     # deriv of inv link
  g <- c(0.35, 0.15, 0.25, 0.2, 0.05, 0.1)      # outcome model parameters
  formula <- "~A*L1 + A*L2"                     # outcome model formula
  ps.formula <- "~L1 + L2"                      # propensity score model formula
  ipw.formula <- "~A"                           # ipw.formula

  # simulate data -----------------------------------------------------------

  set.seed(seed)
  cov.e <- vare
  L1 <- rbinom(n, 1, 0.5)                                  # confounder 1
  L2 <- rnorm(n, 0, 0.16)                                  # confounder 2
  EA <- 0.1 - 0.1*L1 + 0.3*L2                              # E(A|L)
  A <- rnorm(n, EA, sqrt(0.04))                            # exposure
  Astar <- A + rnorm(n, 0, sqrt(cov.e))                    # mismeasured A
  a <- c(0, 1)                                             # grid of exposures
  EY <- inv.link(model.matrix(as.formula(formula)) %*% g)  # mean of outcome
  #hist(EY); range(EY)
  EY[EY < 0] <- 0; EY[EY > 1] <- 1                         # constrain EY
  Y <- rbinom(n, 1, EY)                                    # outcome
  R <- rbinom(n, 1, pi.cc)                                 # c-c sampling
  A[R == 0 & Y == 0] <-
    Astar[R == 0 & Y == 0] <- NA
  dat0 <- data.frame(Y, A, L1, L2, R)                      # oracle data
  datstar <- data.frame(Y, A = Astar, L1, L2, R)           # measured data
  args <- list(formula = formula,                          # arguments for fitting
               ps.formula = ps.formula,
               inv.link = inv.link,
               d.inv.link = d.inv.link)

  # estimate case-cohort weights --------------------------------------------

  pi.cc.hat <- mean(datstar$R[datstar$Y == 0])
  dat0$cc.wts <- datstar$cc.wts <- (1 - datstar$Y) * datstar$R / pi.cc.hat + Y

  # extract estimate and se using delta method ------------------------------

  get.est.se.a <- function(res) {

    est <- diff(tail(unname(res$est), 2))
    vec <- numeric(length(res$est))
    vec[(length(vec)-1):length(vec)] <- c(1, -1)
    se <- sqrt(vec %*% res$var %*% vec)
    bse <- sqrt(vec %*% res$bc.var %*% vec)

    c(est = est,
      se = se,
      bse = bse)
  }

  # estimate dose response curve --------------------------------------------

  args <- list(formula = formula,                       # arguments for fitting
               ps.formula = ps.formula,
               inv.link = inv.link,
               d.inv.link = d.inv.link)

  # naive doubly robust
  dr.naive0 <- fit.dr(data = datstar, args = args, a = a)
  dr.naive <- get.est.se.a(dr.naive0)

  # oracle doubly robust
  dr.oracle0 <- fit.dr(data = dat0, args = args, a = a,
                       start = dr.naive0$est[1:6])
  dr.oracle <- get.est.se.a(dr.oracle0)

  # corrected doubly robust
  dr.mccs0 <- fit.dr.mccs(data = datstar, args = args, a = a,
                          cov.e = cov.e, B = B, mc.seed = mc.seed,
                          start = dr.naive0$est[1:6])
  dr.mccs <- get.est.se.a(dr.mccs0)


  # combine results ---------------------------------------------------------

  res <- c(n, vare, B, seed, pi.cc,
           dr.naive[1], dr.oracle[1], dr.mccs[1],
           dr.naive[2], dr.oracle[2], dr.mccs[2],
           dr.naive[3], dr.oracle[3], dr.mccs[3])

  # return named numeric vector of length: 14
  names(res) <- c(
    "n", "vare", "B", "seed", "pi.cc",
    apply(tidyr::expand_grid(
      c("ghat", "stde", "bste"),
      c("N", "O", "C"),
      1), 1, paste, collapse="."))

  return(res)
}

