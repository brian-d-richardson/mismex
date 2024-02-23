#' Oracle GLM estimating function
#'
#' @param Y a numeric vector, outcome variable
#' @param A a matrix of numbers, possibly multivariate continuous exposure
#'    values with rows corresponding to observations
#' @param L a vector of numbers, confounding variable
#' @param g a numeric vector, coefficients in GLM Y~A+L
#' @param inv.link a function, inverse link function
#' @param d.inv.link a function, derivative of inv.link
#' @param return.sums a logical, indicator for whether a sum of estimating
#'    function values (as opposed to individual values) is to be returned
#'    (default is TRUE)
#'
#' @return individual or summation estimating function values
#'
#' @export
get.psi.glm <- function(data, g,
                        args,
                        return.sums = T) {

  ## unpack arguments
  list2env(args, envir = environment())

  ## extract values from data
  Y <- data$Y
  A <- data[,grepl("A", colnames(data))]
  L <- data[,grepl("L", colnames(data))]

  ## design matrix
  X <- mod.mat(trms = terms(as.formula(formula)),
               data = data)

  ## evaluate estimating equation
  psi <- as.vector((Y - inv.link(X %*% g)) *
                    d.inv.link(X %*% g)) *
    X |>
    `colnames<-`(as.character(1:ncol(X)))

  if (return.sums) {
    return(colSums(psi))
  } else {
    return(psi)
  }
}

#' propensity score estimating function
#'
#' @param A a matrix of numbers, possibly multivariate continuous exposure
#'    values with rows corresponding to observations
#' @param L a vector of numbers, confounding variable
#' @param coef.a.l a vector of numbers, coefficients for the linear model A~L
#' @param var.a.l a vector of numbers, values of the covariance matrix of A|L
#' @param return.sums a logical, indicator for whether a sum of estimating
#'    function values (as opposed to individual values) is to be returned
#'    (default is TRUE)
#'
#' @return individual or summation estimating function values
#'
#' @export
get.psi.ps <- function(data, ps.formula, coef.a.l, var.a.l, return.sums = T) {

  n <- nrow(data)                                   # sample size
  ind.A <- grepl("A", colnames(data))               # exposure columns
  A <- as.matrix(data[,ind.A])
  len.A <- ncol(A)                                  # dimension of exposure
  modmat <- mod.mat(terms(as.formula(ps.formula)),  # PS model matrix
                    data = data)
  len.ps <- ncol(modmat)                            # dimension of PS params

  # residual matrix
  rsd <- A - modmat %*% t(coef.a.l)

  xx <- vapply(1:ncol(modmat), FUN.VAL = rsd,
               FUN = function(x) rsd * modmat[,x])

  psi.ps <- cbind(
    do.call(cbind, lapply(1:len.ps, function(i) rsd * modmat[, i])), # mean                                             # mean
    rsd ^ 2 - matrix(var.a.l, nrow = n, ncol = len.A, byrow = T))    # variance

  if (return.sums) {
    return(colSums((psi.ps)))
  } else {
    return(psi.ps)
  }
}


#' compute standardized IP weights for multivariate normal exposure
#'
#' @inheritParams get.psi.ps
#' @param mean.a a vector of numbers, the marginal mean of the exposure A
#' @param cov.a a matrix of numbers, the marginal covariance of the exposure A
#'
#' @return a vector of standardized IP weights
#'
#' @export
get.SW <- function(data,
                   ps.formula,
                   coef.a.l, var.a.l,
                   mean.a, cov.a) {

  n <- nrow(data) # sample size
  A <- as.matrix(data[,grepl("A", colnames(data))]) # exposure matrix

  # invert covariance of A
  cov.a.inv <- solve(cov.a)

  # set function value type (real or complex)
  fun.val <- ifelse(is.complex(ifelse(is.vector(A), A[1], A[1, 1])),
                    complex(1, 0, 0),
                    0)

  # PS model matrix
  modmat <- mod.mat(terms(as.formula(ps.formula)), data = data)

  # standardized weights
  SW <- sqrt(prod(var.a.l) / det(as.matrix(cov.a))) *
    exp(0.5 *
          vapply(X = 1:n, FUN.VALUE = fun.val, FUN = function(ii) {
            o1 <- A[ii,] - modmat[ii,] %*% t(coef.a.l)
            i1 <- diag(var.a.l ^ -1, nrow = length(var.a.l))
            o2 <- t(A[ii,] - mean.a)
            o1 %*% i1 %*% t(o1) - o2 %*% cov.a.inv %*% t(o2)
          }))

  return(SW)
}


#' Oracle IPW estimating function
#'
#' @inheritParams get.psi.ps
#'
#' @param Y a numeric vector, outcome variable
#' @param g a numeric vector, coefficients in the MSM Y~A
#' @param inv.link a function, link function
#' @param d.inv.link a function, derivative of inv.link
#'
#' @return individual or summation estimating function values
#'
#' @export
get.psi.ipw <- function(data, g, args,
                        mean.a, cov.a,
                        return.sums = T) {

  ## unpack arguments
  list2env(args, envir = environment())

  ## design matrix
  X <- mod.mat(trms = terms(as.formula(formula)),
               data = data)

  ## extract dimensions
  len.msm <- ncol(X)                                         # dim of msm params
  len.ps <- ncol(model.matrix(terms(as.formula(ps.formula)), # PS model params
                              data = data))                  # dim of ps params
  ind.A <- grepl("A", colnames(data))                        # exposure columns
  len.A <- sum(ind.A)

  ## extract PS model params
  coef.a.l <- matrix(g[len.msm + 1:(len.A*len.ps)],
                     ncol = len.ps, byrow = F)
  var.a.l <- exp(tail(g, len.A))

  ## IPW estimating function
  psi <- get.SW(data = data, ps.formula = ps.formula,
                coef.a.l = coef.a.l, var.a.l = var.a.l,
                mean.a = mean.a, cov.a = cov.a) *
    as.vector((data$Y - inv.link(X %*% g[1:len.msm])) *
               d.inv.link(X %*% g[1:len.msm])) *
    X

  if (return.sums) {
    return(colSums(psi))
  } else {
    return(psi)
  }
}

#' create MCCS estimating function
#'
#' @param psi an estimating function
#'
#' @return an MCCS estimating function
#'
#' @export
make.mccs <- function(get.psi, data, args, cov.e, B, mc.seed) {

  ## unpack arguments
  list2env(args, envir = environment())

  ## store dimensions
  n <- nrow(data)                                    # sample size
  ind.A <- grepl("A", colnames(data))                # exposure columns
  len.A <- sum(ind.A)                                # dimension of exposure
  len.est <- ncol(model.matrix(                      # dim of model parameters
    terms(as.formula(formula)), data = data))

  ## MCCS estimating function
  get.psi.mccs <- function(x, return.sums = T) {

    ## set seed before Monte-Carlo sims
    set.seed(mc.seed)

    psi <- vapply(
      X = 1:B,
      FUN = function(b) {

        ## sample n imaginary measurement error values
        e.tilde <- matrix(complex(
          len.A, real = 0,
          imaginary = mvrnorm(n = n,
                              mu = rep(0, len.A),
                              Sigma = cov.e)),
          nrow = n, byrow = F)
        data.tilde <- data
        data.tilde[, ind.A] <- data[, ind.A] + e.tilde

        ## evaluate psi at Astar + e.tilde and take real component
        Re(get.psi(data = data.tilde, g = x,
                   args = args, return.sums = F))
      },
      FUN.VALUE = numeric(n * len.est)) |>
      rowMeans() |>
      matrix(nrow = n, ncol = len.est, byrow = F)

    if (return.sums) {
      return(colSums(psi))
    } else {
      return(psi)
    }
  }
  return(get.psi.mccs)
}
