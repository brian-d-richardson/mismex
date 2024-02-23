#' Fit oracle GLM estimating equation
#'
#' @inheritParams get.psi.glm
#'
#' @param start an optional numeric vector, starting parameter values
#'
#' @return root of GLM estimating function
#'
#' @export
fit.glm <- function(data, args,
                    start = NULL, return.var = TRUE) {

  ## unpack arguments
  list2env(args, envir = environment())

  ## store dimensions
  n <- nrow(data)                               # sample size
  len.est <- ncol(model.matrix(                 # dimension of model parameters
    terms(as.formula(formula)), data = data))

  ## set starting value if not supplied
  if (is.null(start)) {
    start <- rep(0, len.est)
  }

  ## solve oracle IPW equation
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

  return(list(est = root,
              var = evar))
}


#' Fit MCCS GLM estimating equation
#'
#' @param start an optional numeric vector, starting parameter values
#'
#' @return root of MCCS GLM estimating function
#'
#' @export
fit.glm.mccs <- function(data, args,
                         cov.e, B, mc.seed,
                         start = NULL, return.var = TRUE) {

  ## unpack arguments
  list2env(args, envir = environment())

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

  return(list(est = root,
              var = evar))
}
