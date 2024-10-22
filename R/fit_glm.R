#' Fit oracle GLM estimating equation
#'
#' @inheritParams get.psi.glm
#'
#' @param start an optional numeric vector, starting parameter values
#' @param return.var an indicator for whether empirical sandwich variance
#' estimator should be computed, default is TRUE
#' @param return.bcvar an indicator for whether bias-corrected variance
#' estimator should be computed, default is TRUE
#'
#' @return a list of arguments including
#' \itemize{
#' \item{`est`: root of estimating function}
#' \item{`var`: estimated covariance matrix of estimator (if requested)}
#' \item{`bc.var`: bias-corrected estimated covariance matrix of estimator
#' (if requested)}
#' }
#'
#' @export
fit.glm <- function(data, args,
                    start = NULL,
                    return.var = TRUE, return.bcvar = TRUE) {

  ## unpack arguments
  list2env(args, envir = environment())

  ## subset to complete cases if using case-cohort data
  if ("R" %in% colnames(data)) {
    data <- data[data$R == 1 | data$Y == 1,]
  }

  ## add case-cohort weights of 1 if not supplied
  if ( !("cc.wts" %in% colnames(data)) ) {
    data$cc.wts <- 1
  }

  ## store dimensions
  n <- nrow(data)                               # sample size
  len.est <- ncol(model.matrix(                 # dimension of model parameters
    terms(as.formula(formula)), data = data))

  ## set starting value if not supplied
  if (is.null(start)) {
    start <- rep(0, len.est)
  }

  ## solve oracle equation
  root <- tryCatch(
    expr = rootSolve::multiroot(
      f = function(x) get.psi.glm(
        g = x, data = data, args = args),
      start = start)$root,
    warning = function(w) {message(w); rep(NA, len.est)},
    error = function(e) {message(e); rep(NA, len.est)})
  names(root) <- paste0("g.", 0:(len.est - 1))

  # sandwich variance estimate if requested
  evar = matrix(NA, len.est, len.est)
  if (return.var) {
    evar <- tryCatch(
      expr = get.sand.est(
        ghat = root,
        n = n,
        get.psi = function(x) get.psi.glm(
          data = data, g = x, args = args, return.sums = F)),
      warning = function(w) {message(w); matrix(NA, len.est, len.est)},
      error = function(e) {message(e); matrix(NA, len.est, len.est)})
  }

  # bias-corrected sandwich variance estimate if requested
  bc.evar = matrix(NA, len.est, len.est)
  if (return.bcvar) {
    bc.evar <- tryCatch(
      expr = get.sand.est.bc(
        ghat = root,
        n = n,
        get.psi = function(x) get.psi.glm(
          data = data, g = x, args = args, return.sums = F)),
      warning = function(w) {message(w); bc.evar },
      error = function(e) {message(e); bc.evar })
  }

  return(list(est = root,
              var = evar,
              bc.var = bc.evar))
}


#' Fit MCCS GLM estimating equation
#'
#' @inheritParams fit.glm
#' @inheritParams make.mccs
#'
#' @return a list of arguments including
#' \itemize{
#' \item{`est`: root of estimating function}
#' \item{`var`: estimated covariance matrix of estimator (if requested)}
#' \item{`bc.var`: bias-corrected estimated covariance matrix of estimator
#' (if requested)}
#' }
#'
#' @export
fit.glm.mccs <- function(data, args,
                         cov.e, B, mc.seed,
                         start = NULL,
                         return.var = TRUE, return.bcvar = TRUE) {

  ## unpack arguments
  list2env(args, envir = environment())

  ## subset to complete cases if using case-cohort data
  if ("R" %in% colnames(data)) {
    data <- data[data$R == 1 | data$Y == 1,]
  }

  ## add case-cohort weights of 1 if not supplied
  if ( !("cc.wts" %in% colnames(data)) ) {
    data$cc.wts <- 1
  }

  ## store dimensions
  n <- nrow(data)                               # sample size
  len.est <- ncol(model.matrix(                 # dimension of model parameters
    terms(as.formula(formula)), data = data))

  ## set starting value if not supplied
  if (is.null(start)) {
    start <- rep(0, len.est)
  }

  ## get naive estimates to use as starting values in root search
  root.naive <- fit.glm(data = data,
                        args = list(formula = formula,
                                    inv.link = inv.link,
                                    d.inv.link = d.inv.link),
                        start = start, return.var = F)$est

  ## create MCCS GLM estimating function
  get.psi.glm.mccs <- make.mccs(
    get.psi = get.psi.glm, data = data, args = args,
    cov.e = cov.e, B = B, mc.seed = mc.seed)

  ## Solve MCCS GLM equation
  root <- tryCatch(
    expr = rootSolve::multiroot(
      f = get.psi.glm.mccs,
      start = root.naive)$root,
    warning = function(w) {message(w); rep(NA, len.est)},
    error = function(e) {message(e); rep(NA, len.est)})
  names(root) <- paste0("g.", 0:(len.est - 1))

  # sandwich variance estimate if requested
  evar <- matrix(NA, len.est, len.est)
  if (return.var) {
    evar <- tryCatch(
      expr = get.sand.est(
        ghat = root,
        n = n,
        get.psi = function(x) get.psi.glm.mccs(x, return.sums = F)),
      warning = function(w) {message(w); evar},
      error = function(e) {message(e); evar})
  }

  # bias-corrected sandwich variance estimate if requested
  bc.evar = matrix(NA, len.est, len.est)
  if (return.bcvar) {
    bc.evar <- tryCatch(
      expr = get.sand.est.bc(
        ghat = root,
        n = n,
        get.psi = function(x) get.psi.glm.mccs(x, return.sums = F)),
      warning = function(w) {message(w); bc.evar },
      error = function(e) {message(e); bc.evar })
  }

  return(list(est = root,
              var = evar,
              bc.var = bc.evar))
}
