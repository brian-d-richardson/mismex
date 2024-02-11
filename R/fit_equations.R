#' Fit oracle GLM estimating equation
#'
#' @inheritParams get.psi.glm
#'
#' @param start an optional numeric vector, starting parameter values
#'
#' @return root of GLM estimating function
#'
#' @export
fit.glm <- function(Y, A, L, inv.link, d.inv.link,
                    start = NULL, return.var = TRUE) {

  len.a <- ncol(A)              # dimension of A
  len.est <- 2 * (len.a + 1)    # dimension of model parameters
  n <- nrow(A)                  # sample size

  # set starting value if not supplied
  if (is.null(start)) {
    start <- rep(0, len.est)
  }

  # solve oracle IPW equation
  root <- tryCatch(
    expr = rootSolve::multiroot(
      f = function(g) get.psi.glm(g = g, Y = Y, A = A, L = L,
                                  inv.link = inv.link, d.inv.link = d.inv.link),
      start = start)$root,
    warning = function(w) rep(NA, len.est),
    error = function(e) rep(NA, len.est))
  names(root) <- paste0("g.", 0:(len.est - 1))

  # sandwich variance estimate if requested
  evar = matrix(NA, len.est, len.est)
  if (return.var) {
    evar <- tryCatch(
      expr = get.sand.est(
        ghat = root,
        len.a = len.a, n = n,
        get.psi = function(x) get.psi.glm(
          Y = Y, A = A, L = L, g = x,
          inv.link = inv.link, d.inv.link = d.inv.link, return.sums = F)),
      warning = function(w) matrix(NA, len.est, len.est),
      error = function(e) matrix(NA, len.est, len.est))
  }

  return(list(est = root,
              var = evar))
}

#' Fit oracle IPW estimating equation
#'
#' @inheritParams get.psi.ipw
#'
#' @param start an optional numeric vector, starting parameter values
#'
#' @return root of IPW estimating function
#'
#' @export
fit.ipw <- function(Y, A, L,
                    inv.link, d.inv.link,
                    start = NULL,
                    return.var = TRUE,
                    mean.a = NULL, cov.a = NULL,
                    coef.a.l = NULL, var.a.l = NULL) {

  len.a <- ncol(A)            # dimension of A
  len.msm <- len.a + 1        # dimension of MSM parameters
  len.ps <- 3 * len.a         # dimension of propensity score model parameters
  n <- nrow(A)                # sample size

  # compute marginal mean and covariance of A if not supplied
  if (is.null(mean.a)) {mean.a <- colMeans(A)}
  if (is.null(cov.a)) {cov.a <- cov(A)}

  # set starting value if not supplied
  if (is.null(start)) {start <- rep(0, len.msm)}

  # fit propensity score model if not supplied
  if (is.null(coef.a.l) | is.null(var.a.l)) {
    model.a.l <- lm(A ~ L)
    coef.a.l <- t(coef(model.a.l))
    var.a.l <- apply(model.a.l$residuals, 2, var)
  }

  # solve IPW equation
  root <- tryCatch(
    expr = rootSolve::multiroot(
      f = function(x) {
        get.psi.ipw(
          Y = Y, A = A, L = L, g = x[1:(len.a + 1)],
          inv.link = inv.link, d.inv.link = d.inv.link,
          coef.a.l = coef.a.l, var.a.l = var.a.l,
          mean.a = mean.a, cov.a = cov.a) },
      start = start)$root,
    warning = function(w) rep(NA, len.msm),
    error = function(e) rep(NA, len.msm))

  # combine MSM and PS model parameters
  est <- c(root, coef.a.l, log(var.a.l))
  names(est) <- c(
    paste0("g.", 0:len.a),
    paste0("coef.a.l.", 1:(2*len.a)),
    paste0("log.var.a.l", 1:len.a))

  # sandwich variance estimate if requested
  evar <- matrix(NA, len.msm + len.ps, len.msm + len.ps)
  if (return.var) {
    evar <- tryCatch(
      expr = get.sand.est(
        ghat = est,
        len.a = len.a, n = n,
        get.psi = function(x) {
          coef.a.l <- matrix(x[len.a + 1 + 1:(2 * len.a)],
                             nrow = len.a, byrow = F)
          var.a.l <- exp(x[3 * len.a + 1 + 1:len.a])
          cbind(
            get.psi.ipw(
              Y = Y, A = A, L = L, g = x[1:(len.a + 1)],
              inv.link = inv.link, d.inv.link = d.inv.link,
              coef.a.l = coef.a.l, var.a.l = var.a.l,
              mean.a = mean.a, cov.a = cov.a,
              return.sums = F),
            get.psi.ps(
              A = A, L = L,
              coef.a.l = coef.a.l, var.a.l = var.a.l,
              return.sums = F)) }),
      warning = function(w) evar,
      error = function(e) evar)
  }

  return(list(est = est,
              var = evar))
}

#' Fit MCCS GLM estimating equation
#'
#' @inheritParams get.psi.glm.mccs
#'
#' @param start an optional numeric vector, starting parameter values
#'
#' @return root of MCCS GLM estimating function
#'
#' @export
fit.glm.mccs <- function(Y, Astar, L,
                         var.e, B = 10, seed = 123,
                         inv.link, d.inv.link,
                         return.var = TRUE,
                         start = NULL) {

  len.a <- ncol(Astar)          # dimension of A
  len.est <- 2 * (len.a + 1)    # dimension of model parameters
  n <- nrow(Astar)              # sample size

  # set starting value if not supplied
  if (is.null(start)) {
    start <- rep(0, len.est)
  }

  # get naive estimates to use as starting values
  root.naive <- fit.glm(Y = Y, A = Astar, L = L,
                        inv.link = inv.link, d.inv.link = d.inv.link,
                        start = start, return.var = F)$est

  # Solve MCCS GLM equation
  root <- tryCatch(
    expr = rootSolve::multiroot(
      f = function(x) {
        get.psi.glm.mccs(
          Y = Y, Astar = Astar, L = L, g = x,
          inv.link = inv.link, d.inv.link = d.inv.link,
          var.e = var.e, B = B, seed = seed) },
      start = root.naive)$root,
    warning = function(w) rep(NA, len.est),
    error = function(e) rep(NA, len.est))
  names(root) <- paste0("g.", 0:(len.est - 1))

  # sandwich variance estimate if requested
  evar <- matrix(NA, len.est, len.est)
  if (return.var) {
    evar <- tryCatch(
      expr = get.sand.est(
        ghat = root,
        len.a = len.a, n = n,
        get.psi = function(x) get.psi.glm.mccs(
          Y = Y, Astar = Astar, L = L, g = x,
          var.e = var.e, B = B, seed = 123,
          inv.link = inv.link, d.inv.link = d.inv.link, return.sums = F)),
      warning = function(w) evar,
      error = function(e) evar)
  }

  return(list(est = root,
              var = evar))
}

#' Fit MCCS IPW estimating equation
#'
#' @inheritParams get.psi.ipw
#'
#' @param start an optional numeric vector, starting parameter values
#'
#' @return root of MCCS IPW estimating function
#'
#' @export
fit.ipw.mccs <- function(Y, Astar, L,
                         var.e, B = 10, seed = 123,
                         inv.link, d.inv.link,
                         start = NULL,
                         return.var = TRUE,
                         mean.a = NULL, cov.a = NULL,
                         coef.a.l = NULL, var.a.l = NULL) {

  len.a <- ncol(Astar)        # dimension of A
  len.msm <- len.a + 1        # dimension of MSM parameters
  len.ps <- 3 * len.a         # dimension of propensity score model parameters
  n <- nrow(Astar)            # sample size

  # compute marginal mean and covariance of A if not supplied
  if (is.null(mean.a)) {mean.a <- colMeans(Astar)}
  if (is.null(cov.a)) {cov.a <- cov(Astar) - diag(var.e)}

  # set starting value if not supplied
  if (is.null(start)) {start <- rep(0, len.a + 1)}

  # fit propensity score model if not supplied
  if (is.null(coef.a.l)) {
    model.a.l <- lm(Astar ~ L)
    coef.a.l <- t(coef(model.a.l))
    var.a.l <- apply(model.a.l$residuals, 2, var) - var.e
  }

  # get naive estimates to use as starting values
  root.naive <- fit.ipw(Y = Y, A = Astar, L = L,
                        inv.link = inv.link, d.inv.link = d.inv.link,
                        start = start, return.var = F)$est[1:len.msm]

  # Solve MCCS IPW equation
  root <- tryCatch(
    expr =  rootSolve::multiroot(
      f = function(x) {
        get.psi.ipw.mccs(
          Y = Y, Astar = Astar, L = L, g = x,
          var.e = var.e, B = B, seed = seed,
          inv.link = inv.link, d.inv.link = d.inv.link,
          coef.a.l = coef.a.l, var.a.l = var.a.l,
          mean.a = mean.a, cov.a = cov.a) },
      start = root.naive)$root,
    warning = function(w) rep(NA, len.msm),
    error = function(e) rep(NA, len.msm))

  # combine MSM and PS model parameters
  est <- c(root, coef.a.l, log(var.a.l))
  names(est) <- c(
    paste0("g.", 0:len.a),
    paste0("coef.a.l.", 1:(2*len.a)),
    paste0("log.var.a.l", 1:len.a))

  # sandwich variance estimate if requested
  evar <- matrix(NA, len.msm + len.ps, len.msm + len.ps)
  if (return.var) {
    evar <- tryCatch(
      expr = get.sand.est(
        ghat = est,
        len.a = len.a, n = n,
        get.psi = function(x) {
          coef.a.l <- matrix(x[len.a + 1 + 1:(2 * len.a)], nrow = len.a, byrow = F)
          var.a.l <- exp(x[3 * len.a + 1 + 1:len.a])
          cbind(
            get.psi.ipw.mccs(
              Y = Y, Astar = Astar, L = L, g = x[1:((len.a) + 1)],
              var.e = var.e, B = B, seed = 123,
              inv.link = inv.link, d.inv.link = d.inv.link,
              coef.a.l = coef.a.l,
              var.a.l = var.a.l,
              mean.a = mean(Astar), cov.a = cov(Astar) - diag(var.e),
              return.sums = F),
            get.psi.ps(
              A = Astar, L = L,
              coef.a.l = coef.a.l,
              var.a.l = var.a.l + var.e,
              return.sums = F)) }),
      warning = function(w) evar,
      error = function(e) evar)
  }

  return(list(est = est,
              var = evar))
}
















