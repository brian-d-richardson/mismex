#' Fit oracle IPW estimating equation
#'
#' @inheritParams get.psi.ipw
#' @inheritParams get.psi.ps
#' @inheritParams fit.glm
#'
#' @param coef.a.l an optional numeric matrix, coefficients in propensity score
#' model, computed if not specified
#' @param var.a.l an optional numeric vector, variance of A|L in propensity
#' score model, computed if not specified
#' @param mean.a an optional numeric vector, the marginal mean of the exposure
#' A, computed if not specified
#' @param cov.a an optional numeric matrix, the marginal covariance of the
#' exposure A, computed if not specified
#'
#' @return a list of arguments including
#' \itemize{
#' \item{`est`: root of estimating function}
#' \item{`var`: estimated covariance matrix of estimator (if requested)}
#' }
#'
#' @export
fit.ipw <- function(data, args,
                    start = NULL, return.var = TRUE,
                    mean.a = NULL, cov.a = NULL,
                    coef.a.l = NULL, var.a.l = NULL) {

  ## unpack arguments
  list2env(args, envir = environment())

  ## store dimensions
  n <- nrow(data)                                            # sample size
  ind.A <- grepl("A", colnames(data))                        # exposure columns
  A <- as.matrix(data[,ind.A])
  L <- as.matrix(data[,grepl("L", colnames(data))])
  len.A <- sum(ind.A)                                        # dim of exposure
  len.msm <- ncol(model.matrix(                              # dim of msm params
    terms(as.formula(formula)), data = data))
  len.ps <- ncol(model.matrix(                               # dim of ps params
    terms(as.formula(ps.formula)), data = data))

  # compute marginal mean and covariance of A if not supplied
  if (is.null(mean.a)) { mean.a <- colMeans(as.matrix(A)) }
  if (is.null(cov.a)) { if (is.vector(A)) cov.a <- var(A) else cov.a <- cov(A) }

  # set starting value if not supplied
  if (is.null(start)) { start <- rep(0, len.msm) }

  # fit propensity score model
  if (is.null(coef.a.l) | is.null(var.a.l)) {
    model.a.l <- lm(as.formula(paste0("A", ps.formula)), data = data)
    coef.a.l <- t(coef(model.a.l))
    var.a.l <- apply(as.matrix(model.a.l$residuals, ncol = len.a), 2, var)
  }

  # solve IPW equation
  root <- tryCatch(
    expr = rootSolve::multiroot(
      f = function(x) {
        get.psi.ipw(
          data = data, g = c(x, c(coef.a.l), log(var.a.l)),
          args = args, mean.a = mean.a, cov.a = cov.a) },
      start = start)$root,
    warning = function(w) {message(w); rep(NA, len.msm)},
    error = function(e) {message(e); rep(NA, len.msm)})

  # combine MSM and PS model parameters
  est <- c(root, coef.a.l, log(var.a.l))
  names(est) <- c(
    paste0("g.", 0:(len.msm - 1)),
    paste0("coef.a.l.", 1:(len.A*len.ps)),
    paste0("log.var.a.l", 1:len.A))

  # sandwich variance estimate if requested
  evar <- matrix(NA, len.msm + len.ps, len.msm + len.ps)
  if (return.var) {
    evar <- tryCatch(
      expr = get.sand.est(
        ghat = est,
        n = n,
        get.psi = function(x) {
          coef.a.l <- matrix(x[len.msm + 1:(len.A*len.ps)],
                             ncol = len.ps, byrow = F)
          var.a.l <- exp(tail(x, len.A))
          cbind(
            get.psi.ipw(
              data = data, g = x, args = args,
              mean.a = mean.a, cov.a = cov.a,
              return.sums = F),
            get.psi.ps(
              data = data, ps.formula = ps.formula,
              coef.a.l = coef.a.l, var.a.l = var.a.l,
              return.sums = F)) }),
      warning = function(w) {message(w); evar},
      error = function(e) {message(e); evar})
  }

  return(list(est = est,
              var = evar))
}


#' Fit MCCS IPW estimating equation
#'
#' @inheritParams fit.ipw
#' @inheritParams make.mccs
#'
#' @return a list of arguments including
#' \itemize{
#' \item{`est`: root of estimating function}
#' \item{`var`: estimated covariance matrix of estimator (if requested)}
#' }
#'
#' @export
fit.ipw.mccs <- function(data, args,
                         cov.e, B, mc.seed = 123,
                         start = NULL, return.var = TRUE,
                         mean.a = NULL, cov.a = NULL,
                         coef.a.l = NULL, var.a.l = NULL) {

  ## unpack arguments
  list2env(args, envir = environment())

  ## store dimensions
  n <- nrow(data)                                            # sample size
  ind.A <- grepl("A", colnames(data))                        # exposure columns
  A <- as.matrix(data[,ind.A])
  L <- as.matrix(data[,grepl("L", colnames(data))])
  len.A <- sum(ind.A)                                        # dim of exposure
  len.msm <- ncol(model.matrix(                              # dim of msm params
    terms(as.formula(formula)), data = data))
  len.ps <- ncol(model.matrix(terms(as.formula(ps.formula)), # PS model params
                              data = data))                  # dim of ps params
  d.cov.e <- diag(as.matrix(cov.e))                          # cov.e vector

  # set starting value if not supplied
  if (is.null(start)) { start <- rep(0, len.msm) }

  # compute marginal mean and covariance of A if not supplied
  if (is.null(mean.a)) { mean.a <- colMeans(as.matrix(A)) }
  if (is.null(cov.a)) {
    if (is.vector(A)) {
      cov.a <- var(A) - cov.e
    } else {
      cov.a <- cov(A) - cov.e
    }
  }

  # fit propensity score model if not supplied
  if (is.null(coef.a.l)) {
    model.a.l <- lm(as.formula(paste0("A", ps.formula)), data = data)
    coef.a.l <- t(coef(model.a.l))
    var.a.l <- apply(as.matrix(model.a.l$residuals, ncol = len.a), 2, var) -
      d.cov.e
  }

  ## get naive estimates to use as starting values
  root.naive <- fit.ipw(data = data, args = args,
                        start = start, return.var = F)$est[1:len.msm]

  ## create MCCS IPW estimating function
  get.psi.ipw.mccs <- make.mccs(
    get.psi = function(data, g, args, return.sums = T) {
      get.psi.ipw(data = data, args = args,
                  g = g,
                  mean.a = mean.a, cov.a = cov.a,
                  return.sums = return.sums) },
    data = data, args = args,
    cov.e = cov.e, B = B, mc.seed = mc.seed)

  # Solve MCCS IPW equation
  root <- tryCatch(
    expr = rootSolve::multiroot(
      f = function(xx) get.psi.ipw.mccs(x = c(xx, c(coef.a.l), log(var.a.l))),
      start = root.naive)$root,
    warning = function(w) {message(w); rep(NA, len.msm)},
    error = function(e) {message(e); rep(NA, len.msm)})

  # combine MSM and PS model parameters
  est <- c(root, coef.a.l, log(var.a.l))
  names(est) <- c(
    paste0("g.", 0:(len.msm - 1)),
    paste0("coef.a.l.", 1:(len.A*len.ps)),
    paste0("log.var.a.l", 1:len.A))

  # sandwich variance estimate if requested
  evar <- matrix(NA, len.msm + len.ps, len.msm + len.ps)
  if (return.var) {
    evar <- tryCatch(
      expr = get.sand.est(
        ghat = est,
        n = n,
        get.psi = function(x) {
          cbind(
            get.psi.ipw.mccs(x = x, return.sums = F),
            get.psi.ps(
              data = data, ps.formula = ps.formula,
              coef.a.l = matrix(x[len.msm + 1:(len.A*len.ps)],
                                ncol = len.ps, byrow = F),
              var.a.l = exp(tail(x, len.A)),# + d.cov.e,
              return.sums = F)) }),
      warning = function(w) {message(w); evar},
      error = function(e) {message(e); evar})
  }

  return(list(est = est,
              var = evar))
}

