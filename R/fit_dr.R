#' Fit oracle doubly robust estimating equation
#'
#' @inheritParams get.psi.glm
#' @inheritParams get.psi.ipw
#' @inheritParams fit.glm
#' @param args a list of arguments including
#' \itemize{
#' \item{`inv.link`: a function, inverse link function}
#' \item{`d.inv.link`: a function, derivative of inv.link}
#' \item{`formula`: a character string of outcome model formula}
#' \item{`ps.formula`: a character string of propensity model formula}
#' }
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
fit.dr <- function(data, args, a,
                   start = NULL,
                   return.var = TRUE,
                   return.bcvar = TRUE,
                   mean.a = NULL, cov.a = NULL,
                   coef.a.l = NULL, var.a.l = NULL) {

  #start = NULL; return.var = TRUE; mean.a = NULL; cov.a = NULL; coef.a.l = NULL; var.a.l = NULL

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
  len.a <- nrow(as.matrix(a))                   # size of exposure grid
  ind.A <- grepl("A", colnames(data))
  A <- as.matrix(data[,ind.A])
  len.A <- sum(ind.A)                           # dimension of exposure
  len.ps <- ncol(model.matrix(                  # dimension of ps parameters
    terms(as.formula(ps.formula)), data = data))
  L <- as.matrix(data[,grepl("L", colnames(data))])

  # set starting value if not supplied
  if (is.null(start)) { start <- rep(0, len.est) }

  # compute marginal mean and covariance of A if not supplied
  if (is.null(mean.a)) { mean.a <- colMeans(as.matrix(A)) }
  if (is.null(cov.a)) { if (is.vector(A)) cov.a <- var(A) else cov.a <- cov(A) }

  # fit propensity score model if not supplied
  if (is.null(coef.a.l) | is.null(var.a.l)) {
    model.a.l <- lm(as.formula(paste0("A", ps.formula)),
                    data = data, weights = cc.wts)
    coef.a.l <- t(coef(model.a.l))
    var.a.l <- apply(as.matrix(model.a.l$residuals, ncol = len.a), 2, var)
  }

  # fit weighted outcome model
  root <- fit.ipw(data = data, args = args, start = start,
                  return.var = F, return.bcvar = F,
                  coef.a.l = coef.a.l, var.a.l = var.a.l)$est
  outcome.params <- head(root, len.est)

  # estimate E{Y(a)} for each supplied a value
  EYa <- vapply(X = 1:len.a,
                FUN.VALUE = 0,
                FUN = function(ai) {
                  if (is.vector(a)) {
                    aa <- a[ai]
                    a.mod.mat <- mod.mat(
                      terms(as.formula(formula)),
                      data = data.frame(
                        A = do.call("rbind", replicate(n, aa, simplify = F)),
                        L))
                  } else {
                    aa <- a[ai,]
                    a.mod.mat <- mod.mat(
                      terms(as.formula(formula)),
                      data = data.frame(
                        do.call("rbind", replicate(n, aa, simplify = F)),
                        L))
                  }
                  mean(inv.link(a.mod.mat %*% outcome.params))
                })
  names(EYa) <- paste0("EYa.", 1:len.a)
  ghat <- c(root, EYa)

  # sandwich variance estimates including PS model if requested
  evar <- as.data.frame(matrix(NA, length(ghat), length(ghat)))
  if (return.var) {
    evar <- tryCatch(
      expr = get.sand.est(
        ghat = ghat,
        n = n,
        get.psi = function(x) {
          ght <- head(x, -len.a)                                  # IPW params
          ght.out <- head(x, len.est)                             # Y|A,L coefs
          coef.a.l <- matrix(                                     # A|L coefs
            x[len.est + 1:(len.A*len.ps)],
            ncol = len.ps, byrow = F)
          var.a.l <- exp(x[len.est + (len.A*len.ps) + 1:len.A])   # Var(A|L)
          mean.a <- x[len.est + (len.A*len.ps) + len.A + 1:len.A] # E(A)
          if (len.A == 1) {
            cov.a <- x[len.est + (len.A*len.ps) + 2*len.A + 1]
          } else {
            cov.a <- matrix(0, len.A, len.A)
            cov.a[upper.tri(cov.a, diag = T)] <-                    # Cov(A)
              x[len.est + (len.A*len.ps) + 2*len.A +
                  1:(len.A * (len.A + 1) / 2)]
            cov.a <- cov.a + t(cov.a) - diag(diag(cov.a))
          }
          EYa <- tail(x, len.a)
          cbind(
            # weighted outcome model
            get.psi.ipw(data = data, g = ght, args = args,
                        mean.a = mean.a, cov.a = cov.a,
                        return.sums = F),
            # PS weights
            get.psi.ps(
              data = data, ps.formula = ps.formula,
              coef.a.l = coef.a.l, var.a.l = var.a.l,
              return.sums = F),
            # PS numerator
            get.psi.ps.num(
              data = data,
              cov.a = cov.a, mean.a = mean.a,
              return.sums = F),
            # E{Y(a)}
            vapply(X = 1:len.a,
                   FUN.VALUE = numeric(n),
                   FUN = function(ai) {
                     if (is.vector(a)) {
                       aa <- a[ai]
                       a.mod.mat <- mod.mat(
                         terms(as.formula(formula)),
                         data = data.frame(
                           A = do.call("rbind", replicate(n, aa, simplify = F)),
                           L))
                     } else {
                       aa <- a[ai,]
                       a.mod.mat <- mod.mat(
                         terms(as.formula(formula)),
                         data = data.frame(
                           do.call("rbind", replicate(n, aa, simplify = F)),
                           L))
                     }
                     EYa[ai] - inv.link(a.mod.mat %*% ght.out)
                   }))}),
      warning = function(w) {message(w); evar },
      error = function(e) {message(e); evar })
  }
  colnames(evar) <- names(ghat)

  # bias-corrected variance if requested
  bc.evar = as.data.frame(matrix(NA, length(ghat), length(ghat)))
  if (return.bcvar) {
    bc.evar <- tryCatch(
      expr = get.sand.est.bc(
        ghat = ghat,
        n = n,
        get.psi = function(x) {
          ght <- head(x, -len.a)                                  # IPW params
          ght.out <- head(x, len.est)                             # Y|A,L coefs
          coef.a.l <- matrix(                                     # A|L coefs
            x[len.est + 1:(len.A*len.ps)],
            ncol = len.ps, byrow = F)
          var.a.l <- exp(x[len.est + (len.A*len.ps) + 1:len.A])   # Var(A|L)
          mean.a <- x[len.est + (len.A*len.ps) + len.A + 1:len.A] # E(A)
          if (len.A == 1) {
            cov.a <- x[len.est + (len.A*len.ps) + 2*len.A + 1]
          } else {
            cov.a <- matrix(0, len.A, len.A)
            cov.a[upper.tri(cov.a, diag = T)] <-                    # Cov(A)
              x[len.est + (len.A*len.ps) + 2*len.A +
                  1:(len.A * (len.A + 1) / 2)]
            cov.a <- cov.a + t(cov.a) - diag(diag(cov.a))
          }
          EYa <- tail(x, len.a)
          cbind(
            # weighted outcome model
            get.psi.ipw(data = data, g = ght, args = args,
                        mean.a = mean.a, cov.a = cov.a,
                        return.sums = F),
            # PS weights
            get.psi.ps(
              data = data, ps.formula = ps.formula,
              coef.a.l = coef.a.l, var.a.l = var.a.l,
              return.sums = F),
            # PS numerator
            get.psi.ps.num(
              data = data,
              cov.a = cov.a, mean.a = mean.a,
              return.sums = F),
            # E{Y(a)}
            vapply(X = 1:len.a,
                   FUN.VALUE = numeric(n),
                   FUN = function(ai) {
                     if (is.vector(a)) {
                       aa <- a[ai]
                       a.mod.mat <- mod.mat(
                         terms(as.formula(formula)),
                         data = data.frame(
                           A = do.call("rbind", replicate(n, aa, simplify = F)),
                           L))
                     } else {
                       aa <- a[ai,]
                       a.mod.mat <- mod.mat(
                         terms(as.formula(formula)),
                         data = data.frame(
                           do.call("rbind", replicate(n, aa, simplify = F)),
                           L))
                     }
                     EYa[ai] - inv.link(a.mod.mat %*% ght.out)
                   }))}),
      warning = function(w) {message(w); bc.evar },
      error = function(e) {message(e); bc.evar })
  }
  colnames(bc.evar) <- names(ghat)

  return(list(est = ghat,
              var = evar,
              bc.var = bc.evar))
}


#' Fit MCCS doubly robust estimating equation
#'
#' @inheritParams fit.dr
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
fit.dr.mccs <- function(data, args, a,
                        cov.e, B, mc.seed = 123,
                        start = NULL,
                        return.var = TRUE,
                        return.bcvar = TRUE,
                        mean.a = NULL, cov.a = NULL,
                        coef.a.l = NULL, var.a.l = NULL) {

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
  len.a <- nrow(as.matrix(a))                   # size of exposure grid
  ind.A <- grepl("A", colnames(data))
  A <- data[,ind.A]                             # exposure
  len.A <- sum(ind.A)                           # dimension of exposure
  len.ps <- ncol(model.matrix(                  # dimension of ps parameters
    terms(as.formula(ps.formula)), data = data))
  L <- data[, grepl("L", colnames(data))]       # covariates
  d.cov.e <- diag(as.matrix(cov.e))             # cov.e vector

  # set starting value if not supplied
  if (is.null(start)) { start <- rep(0, len.est) }

  # compute marginal mean and covariance of A if not supplied
  if (is.null(mean.a)) { mean.a <- colMeans(as.matrix(A)) }
  if (is.null(cov.a)) {
    if (is.vector(A)) {
      cov.a <- var(A) - cov.e
    } else {
      cov.a <- cov(A) - cov.e
    }
  }

  # fit weighted outcome model
  root <- fit.ipw.mccs(data = data, args = args, start = start,
                       return.var = F, return.bcvar = F,
                       cov.e = cov.e, B = B, mc.seed = mc.seed,
                       mean.a = mean.a, cov.a = cov.a,
                       coef.a.l = coef.a.l, var.a.l = var.a.l)$est
  outcome.params <- head(root, len.est)

  # estimate E{Y(a)} for each supplied a value
  EYa <- vapply(X = 1:len.a,
                FUN.VALUE = 0,
                FUN = function(ai) {
                  if (is.vector(a)) {
                    aa <- a[ai]
                    a.mod.mat <- mod.mat(
                      terms(as.formula(formula)),
                      data = data.frame(
                        A = do.call("rbind", replicate(n, aa, simplify = F)),
                        L))
                  } else {
                    aa <- a[ai,]
                    a.mod.mat <- mod.mat(
                      terms(as.formula(formula)),
                      data = data.frame(
                        do.call("rbind", replicate(n, aa, simplify = F)),
                        L))
                    }
                  mean(inv.link(a.mod.mat %*% outcome.params))
                })
  names(EYa) <- paste0("EYa.", 1:len.a)
  ghat <- c(root, EYa)

  ## create MCCS IPW estimating function
  get.psi.ipw.mccs <- make.mccs(
    get.psi = function(data, g, args, return.sums = T) {
      get.psi.ipw(data = data, args = args, g = g,
                  mean.a = mean.a, cov.a = cov.a,
                  return.sums = return.sums) },
    data = data, args = args,
    cov.e = cov.e, B = B, mc.seed = mc.seed)

  # sandwich variance estimates including PS model if requested
  evar = as.data.frame(matrix(NA, length(ghat), length(ghat)))
  if (return.var) {
    evar <- tryCatch(
      expr = get.sand.est(
        ghat = ghat,
        n = n,
        get.psi = function(x) {
          ght <- head(x, -len.a)                                  # IPW params
          ght.out <- head(x, len.est)                             # Y|A,L coefs
          coef.a.l <- matrix(                                     # A|L coefs
            x[len.est + 1:(len.A*len.ps)],
            ncol = len.ps, byrow = F)
          var.a.l <- exp(x[len.est + (len.A*len.ps) + 1:len.A])   # Var(A|L)
          + d.cov.e
          mean.a <- x[len.est + (len.A*len.ps) + len.A + 1:len.A] # E(A)
          if (len.A == 1) {
            cov.a <- x[len.est + (len.A*len.ps) + 2*len.A + 1]
          } else {
            cov.a <- matrix(0, len.A, len.A)
            cov.a[upper.tri(cov.a, diag = T)] <-                  # Cov(A)
              x[len.est + (len.A*len.ps) + 2*len.A +
                  1:(len.A * (len.A + 1) / 2)]
            cov.a <- cov.a + t(cov.a) - diag(diag(cov.a))
          }
          cov.a <- cov.a + cov.e
          EYa <- tail(x, len.a)
          cbind(
            # weighted outcome model
            get.psi.ipw.mccs(x = ght, return.sums = F),
            # PS weights
            get.psi.ps(
              data = data, ps.formula = ps.formula,
              coef.a.l = coef.a.l, var.a.l = var.a.l,
              return.sums = F),
            # PS numerator
            get.psi.ps.num(
              data = data,
              cov.a = cov.a, mean.a = mean.a,
              return.sums = F),
            # E{Y(a)}
            vapply(X = 1:len.a,
                   FUN.VALUE = numeric(n),
                   FUN = function(ai) {
                     if (is.vector(a)) {
                       aa <- a[ai]
                       a.mod.mat <- mod.mat(
                         terms(as.formula(formula)),
                         data = data.frame(
                           A = do.call("rbind", replicate(n, aa, simplify = F)),
                           L))
                     } else {
                       aa <- a[ai,]
                       a.mod.mat <- mod.mat(
                         terms(as.formula(formula)),
                         data = data.frame(
                           do.call("rbind", replicate(n, aa, simplify = F)),
                           L))
                     }
                     EYa[ai] - inv.link(a.mod.mat %*% ght.out)
                   }))}),
      warning = function(w) {message(w); evar },
      error = function(e) {message(e); evar })
  }
  colnames(evar) <- names(ghat)

  # bias-corrected variance if requested
  bc.evar = as.data.frame(matrix(NA, length(ghat), length(ghat)))
  if (return.bcvar) {
    bc.evar <- tryCatch(
      expr = get.sand.est.bc(
        ghat = ghat,
        n = n,
        get.psi = function(x) {
          ght <- head(x, -len.a)                                  # IPW params
          ght.out <- head(x, len.est)                             # Y|A,L coefs
          coef.a.l <- matrix(                                     # A|L coefs
            x[len.est + 1:(len.A*len.ps)],
            ncol = len.ps, byrow = F)
          var.a.l <- exp(x[len.est + (len.A*len.ps) + 1:len.A])   # Var(A|L)
          + d.cov.e
          mean.a <- x[len.est + (len.A*len.ps) + len.A + 1:len.A] # E(A)
          if (len.A == 1) {
            cov.a <- x[len.est + (len.A*len.ps) + 2*len.A + 1]
          } else {
            cov.a <- matrix(0, len.A, len.A)
            cov.a[upper.tri(cov.a, diag = T)] <-                  # Cov(A)
              x[len.est + (len.A*len.ps) + 2*len.A +
                  1:(len.A * (len.A + 1) / 2)]
            cov.a <- cov.a + t(cov.a) - diag(diag(cov.a))
          }
          cov.a <- cov.a + cov.e
          EYa <- tail(x, len.a)
          cbind(
            # weighted outcome model
            get.psi.ipw.mccs(x = ght, return.sums = F),
            # PS weights
            get.psi.ps(
              data = data, ps.formula = ps.formula,
              coef.a.l = coef.a.l, var.a.l = var.a.l,
              return.sums = F),
            # PS numerator
            get.psi.ps.num(
              data = data,
              cov.a = cov.a, mean.a = mean.a,
              return.sums = F),
            # E{Y(a)}
            vapply(X = 1:len.a,
                   FUN.VALUE = numeric(n),
                   FUN = function(ai) {
                     if (is.vector(a)) {
                       aa <- a[ai]
                       a.mod.mat <- mod.mat(
                         terms(as.formula(formula)),
                         data = data.frame(
                           A = do.call("rbind", replicate(n, aa, simplify = F)),
                           L))
                     } else {
                       aa <- a[ai,]
                       a.mod.mat <- mod.mat(
                         terms(as.formula(formula)),
                         data = data.frame(
                           do.call("rbind", replicate(n, aa, simplify = F)),
                           L))
                     }
                     EYa[ai] - inv.link(a.mod.mat %*% ght.out)
                   }))}),
      warning = function(w) {message(w); bc.evar },
      error = function(e) {message(e); bc.evar })
  }
  colnames(bc.evar) <- names(ghat)


  return(list(est = ghat,
              var = evar,
              bc.var = bc.evar))
}


