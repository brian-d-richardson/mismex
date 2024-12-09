#' Simulation 4: DR with model misspecification
#'
#' @inheritParams sim1.gfmla.nonlinear.coarse
#'
#' @return a data frame with the following columns
#' \itemize{
#' \item{n}
#' \item{B}
#' \item{vare}
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
sim4.dr.misspec <- function(
    n,
    vare,
    B,
    seed) {

  # for troubleshooting -----------------------------------------------------

  #library(MASS); library(devtools); load_all();
  #n = 2000; vare = 0.02; B = 80; seed = 1;

  # define parameters -------------------------------------------------------

  mc.seed <- 123                                # MCCS seed
  inv.link <- inv.ident                         # inverse link
  d.inv.link <- d.inv.ident                     # deriv of inv link
  g <- c(0.35, 0.15, 0.25, 0.2, 0.05, 0.1)      # outcome model parameters
  formula <- "~A*L1 + A*L2"                     # outcome model formula
  ps.formula <- "~L1 + L2"                      # propensity score model formula
  ipw.formula <- "~A"                           # ipw.formula
  formula.inc <- "~A*L2"                        # incorrect outcome model
  ps.formula.inc <- "~L2"                       # incorrect propensity score model

  # simulate data -----------------------------------------------------------

  set.seed(seed)
  cov.e <- vare
  L1 <- rbinom(n, 1, 0.5)                                  # confounder 1
  L2 <- rnorm(n, 0, 0.16)                                  # confounder 2
  EA <- 0.1 - 0.1*L1 + 0.3*L2                              # E(A|L)
  A <- rnorm(n, EA, sqrt(0.04))                            # exposure
  Astar <- A + rnorm(n, 0, sqrt(cov.e))                     # mismeasured A
  a <- c(0, 1)                                             # grid of exposures
  EY <- inv.link(model.matrix(as.formula(formula)) %*% g)  # mean of outcome
  #hist(EY); range(EY)
  EY[EY < 0] <- 0; EY[EY > 1] <- 1                         # constrain EY
  Y <- rbinom(n, 1, EY)                                    # outcome
  dat0 <- data.frame(Y, A, L1, L2)                         # oracle data
  datstar <- data.frame(Y, A = Astar, L1, L2)                  # measured data
  args <- list(formula = formula,                          # arguments for fitting
               ps.formula = ps.formula,
               inv.link = inv.link,
               d.inv.link = d.inv.link)

  #plot(L1, A); plot(L2, A); plot(A, Y);
  #var(A) / var(Astar)

  #library(ggplot2)
  #dat0 %>%
  #  mutate(L2grp = factor(ntile(L2, 4)),
  #         L1 = factor(L1)) %>%
  #  ggplot(aes(fill = L2grp,
  #             y = A)) +
  #  geom_boxplot() +
  #  facet_grid(~L1)

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

  # compute estimates and std errors of MSM coefficient for a ---------------

  get.est.se.a <- function(res, res.list) {
    name <- strsplit(res, "[.]")[[1]]

    # for IPW, extract estimate and std error for coefficient of a
    if (name[1] == "ipw") {
      est <- unname(res.list[[res]]$est[2])
      se <- sqrt(diag(res.list[[res]]$var)[2])
      bse <- sqrt(diag(res.list[[res]]$bc.var)[2])

      # for g-fmla and double robust, use delta method on E{Y(1)}, E{Y(0)}
    } else {
      est <- diff(tail(unname(res.list[[res]]$est), 2))
      vec <- numeric(length(res.list[[res]]$est))
      vec[(length(vec)-1):length(vec)] <- c(1, -1)
      se <- sqrt(vec %*% res.list[[res]]$var %*% vec)
      bse <- sqrt(vec %*% res.list[[res]]$bc.var %*% vec)
    }
    c(method = name[1],
      type = name[2],
      est = est,
      se = se,
      bse = bse)
  }

  # (naive, oracle, rc, mccs) x (gfmla, ipw, dr) estimates given model specs
  est.all <- function(ps.formula, formula) {

    # store results in list
    res.list = list()

    # length of outcome model params
    len.out <- ncol(model.matrix(as.formula(formula), data = dat0))

    # g-formula
    gfmla.args <- list(formula = formula,
                       inv.link = inv.link,
                       d.inv.link = d.inv.link)
    res.list[["gfmla.naive"]] <- fit.gfmla(
      data = datstar, args = gfmla.args, a = c(0, 1))
    res.list[["gfmla.oracle"]] <- fit.gfmla(
      data = dat0, args = gfmla.args, a = c(0, 1),
      start = res.list[["gfmla.naive"]]$est[1:len.out])
    res.list[["gfmla.rc"]] <- fit.gfmla(
      data = datrc, args = gfmla.args, a = c(0, 1),
      start = res.list[["gfmla.naive"]]$est[1:len.out])
    res.list[["gfmla.mccs"]] <- fit.gfmla.mccs(
      data = datstar, args = gfmla.args, a = c(0, 1),
      cov.e = cov.e, B = B, mc.seed = mc.seed,
      start = res.list[["gfmla.naive"]]$est[1:len.out])

    # ipw
    ipw.args <- list(formula = ipw.formula,
                     ps.formula = ps.formula,
                     inv.link = inv.link,
                     d.inv.link = d.inv.link)
    res.list[["ipw.naive"]] <- fit.ipw(
      data = datstar, args = ipw.args,
      start = res.list[["gfmla.naive"]]$est[1:2])
    res.list[["ipw.oracle"]] <- fit.ipw(
      data = dat0, args = ipw.args,
      start = res.list[["ipw.naive"]]$est[1:2])
    res.list[["ipw.rc"]] <- fit.ipw(
      data = datrc, args = ipw.args,
      start = res.list[["ipw.naive"]]$est[1:2])
    res.list[["ipw.mccs"]] <- fit.ipw.mccs(
      data = datstar, args = ipw.args,
      cov.e = cov.e, B = B, mc.seed = mc.seed,
      start = res.list[["ipw.naive"]]$est[1:2])

    # dr
    dr.args <- list(formula = formula,
                    ps.formula = ps.formula,
                    inv.link = inv.link,
                    d.inv.link = d.inv.link)
    res.list[["dr.naive"]] <- fit.dr(
      data = datstar, args = dr.args, a = c(0, 1),
      start = res.list[["gfmla.naive"]]$est[1:len.out])
    res.list[["dr.oracle"]] <- fit.dr(
      data = dat0, args = dr.args, a = c(0, 1),
      start = res.list[["dr.naive"]]$est[1:len.out])
    res.list[["dr.rc"]] <- fit.dr(
      data = datrc, args = dr.args, a = c(0, 1),
      start = res.list[["dr.naive"]]$est[1:len.out])
    res.list[["dr.mccs"]] <- fit.dr.mccs(
      data = datstar, args = dr.args, a = c(0, 1),
      cov.e = cov.e, B = B, mc.seed = mc.seed,
      start = res.list[["dr.naive"]]$est[1:len.out])

    dat <- as.data.frame(t(vapply(
      X = names(res.list),
      FUN = function(res) get.est.se.a(res = res, res.list = res.list),
      FUN.VALUE = character(5)))) |>
      mutate_at(c("est", "se", "bse"), as.numeric)

    return(dat)
  }


  # combine results ---------------------------------------------------------

  # (00) both models correct
  res.00 <- est.all(ps.formula = ps.formula,
                    formula = formula)

  # (10) PS incorrect, outcome correct
  res.10 <- est.all(ps.formula = ps.formula.inc,
                    formula = formula)

  # (01) PS correct, outcome incorrect
  res.01 <- est.all(ps.formula = ps.formula,
                    formula = formula.inc)

  # (11) both models incorrect
  res.11 <- est.all(ps.formula = ps.formula.inc,
                    formula = formula.inc)

  # combine results (36 x 11 data frame)
  res <- cbind(n = n, B = B, vare = vare, seed = seed,
               rbind(cbind(ps = 0, out = 0, res.00),
                     cbind(ps = 1, out = 0, res.10),
                     cbind(ps = 0, out = 1, res.01),
                     cbind(ps = 1, out = 1, res.11)))

  return(res)
}

