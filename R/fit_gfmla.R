#' Fit oracle G-formula estimating equation
#'
#' @inheritParams fit.glm
#'
#' @param a exposure value at which to estimate mean potential outcome
#'
#' @return a list of arguments including
#' \itemize{
#' \item{`est`: root of estimating function}
#' \item{`var`: estimated covariance matrix of estimator (if requested)}
#' }
#'
#' @export
fit.gfmla <- function(data, args, a,
                      start = NULL, return.var = TRUE) {

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

  ## store values
  n <- nrow(data)                               # sample size
  len.est <- ncol(model.matrix(                 # dimension of model parameters
    terms(as.formula(formula)), data = data))
  len.a <- length(a)                            # number of exposure values
  L <- data[, grepl("L", colnames(data))]       # covariates

  # fit outcome model
  root <- fit.glm(data = data, args = args, return.var = F)$est

  # estimate E{Y(a)} for each supplied a value
  EYa <- vapply(X = a,
         FUN.VALUE = 0,
         FUN = function(aa) {
           a.mod.mat <- mod.mat(terms(as.formula(formula)),
                                data = data.frame(A = aa, L))
           mean(inv.link(a.mod.mat %*% root))
         })
  names(EYa) <- paste0("EYa.", 1:len.a)
  ghat = c(root, EYa)

  # sandwich variance estimate if requested
  evar = matrix(NA, len.est, len.est)
  if (return.var) {
    evar <- tryCatch(
      expr = get.sand.est(
        ghat = ghat,
        n = n,
        get.psi = function(x) {
          ght <- head(x, -len.a)
          EYa <- tail(x, len.a)
          cbind(
            get.psi.glm(
              data = data, g = ght, args = args, return.sums = F),
            vapply(X = 1:len.a,
                   FUN.VALUE = numeric(n),
                   FUN = function(aa) {
                     a.mod.mat <- mod.mat(terms(as.formula(formula)),
                                          data = data.frame(A = a[aa], L))
                     EYa[aa] - inv.link(a.mod.mat %*% ght)
                   }))}),
      warning = function(w) {message(w); matrix(NA, len.est, len.est)},
      error = function(e) {message(e); matrix(NA, len.est, len.est)})
  }

  return(list(est = ghat,
              var = evar))
}


#' Fit MCCS G-formula estimating equation
#'
#' @inheritParams fit.gfmla
#' @inheritParams make.mccs
#'
#' @return a list of arguments including
#' \itemize{
#' \item{`est`: root of estimating function}
#' \item{`var`: estimated covariance matrix of estimator (if requested)}
#' }
#'
#' @export
fit.gfmla.mccs <- function(data, args, a,
                           cov.e, B, mc.seed,
                           start = NULL,
                           return.var = TRUE) {

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

  ## store values
  data <- data[data$R == 1 | data$Y == 1,]     # complete data
  n <- nrow(dat)                               # sample size
  len.est <- ncol(model.matrix(                # dimension of model parameters
    terms(as.formula(formula)), data = data))
  len.a <- length(a)                            # number of exposure values
  L <- data[, grepl("L", colnames(data))]      # covariates

  ## fit outcome model
  root <- fit.glm.mccs(data = data, args = args,
                       cov.e = cov.e, B = B, mc.seed = mc.seed,
                       return.var = F)$est

  ## estimate E{Y(a)} for each supplied a value
  EYa <- vapply(X = a,
                FUN.VALUE = 0,
                FUN = function(aa) {
                  a.mod.mat <- mod.mat(terms(as.formula(formula)),
                                       data = data.frame(A = aa, L))
                  mean(inv.link(a.mod.mat %*% root))
                })
  names(EYa) <- paste0("EYa.", 1:len.a)
  ghat = c(root, EYa)

  ## create MCCS GLM estimating function
  get.psi.glm.mccs <- make.mccs(
    get.psi = get.psi.glm, data = data, args = args,
    cov.e = cov.e, B = B, mc.seed = mc.seed)

  # sandwich variance estimate if requested
  evar = matrix(NA, len.est + len.a, len.est + len.a)
  if (return.var) {
    evar <- tryCatch(
      expr = get.sand.est(
        ghat = ghat,
        n = n,
        get.psi = function(x) {
          ght <- head(x, -length(a))
          EYa <- tail(x, length(a))
          cbind(
            get.psi.glm.mccs(x = ght, return.sums = F),
            vapply(X = 1:length(a),
                   FUN.VALUE = numeric(n),
                   FUN = function(aa) {
                     a.mod.mat <- mod.mat(terms(as.formula(formula)),
                                          data = data.frame(A = a[aa], L))
                     EYa[aa] - inv.link(a.mod.mat %*% ght)
                   }))}),
      warning = function(w) {message(w); evar},
      error = function(e) {message(e); evar})
  }

  return(list(est = ghat,
              var = evar))
}

