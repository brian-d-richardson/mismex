#' propensity score estimating function
#'
#' @param A ...
#'
#' @return individual or summation estimating function values
#'
#' @export
get.psi.ps <- function(A, L, coef.a.l, var.a.l, return.sums = T) {

  psi.ps <- cbind(
    A - cbind(1, L) %*% t(coef.a.l),                       # mean
    (A - cbind(1, L) %*% t(coef.a.l)) * L,
    (A - cbind(1, L) %*% t(coef.a.l)) ^ 2 -                # variance
      matrix(var.a.l, nrow = nrow(A), ncol = ncol(A), byrow = T))

  if (return.sums) {
    return(colSums((psi.ps)))
  } else {
    return(psi.ps)
  }
}


#' compute standardized weights for (multivariate) normal exposure
#'
#' @param A ...
#'
#' @return individual or summation estimating function values
#'
#' @export
get.SW <- function(A, L,
                   coef.a.l, var.a.l,
                   mean.a, cov.a) {

  # invert covariance of A
  cov.a.inv <- solve(cov.a)

  # set function value type (real or complex)
  fun.val <- ifelse(is.complex(A),
                    complex(1, 0, 0),
                    0)

  SW <- sqrt(prod(var.a.l) / det(cov.a)) *
    vapply(X = 1:nrow(A),
           FUN = function(ii) {
             exp(0.5 * (

               (A[ii,] - cbind(1, L[ii]) %*% t(coef.a.l)) %*%
                 diag(var.a.l ^ -1) %*%
                 t((A[ii,] - cbind(1, L[ii]) %*% t(coef.a.l))) -

                 t(A[ii,] - mean.a) %*% cov.a.inv %*% (A[ii,] - mean.a))) },
           FUN.VALUE = fun.val)

  return(SW)
}


#' Oracle GLM estimating function
#'
#' @param Y outcome, a numeric vector
#' @param A ...
#'
#' @return individual or summation estimating function values
#'
#' @export
get.psi.glm <- function(Y, X, g, inv.link, d.inv.link,
                        return.sums = T) {

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


#' Oracle IPW estimating function
#'
#' @param Y outcome, a numeric vector
#' @param A ...
#'
#' @return individual or summation estimating function values
#'
#' @export
get.psi.ipw <- function(Y, A, L, g, inv.link, d.inv.link,
                        coef.a.l, var.a.l, mean.a, cov.a,
                        return.sums = T) {

  psi <- get.SW(A = A, L = L,
                coef.a.l = coef.a.l, var.a.l = var.a.l,
                mean.a = mean.a, cov.a = cov.a) *
    as.vector((Y - inv.link(cbind(1, A) %*% g)) *
               d.inv.link(cbind(1, A) %*% g)) *
    cbind(1, A) |>
    `colnames<-`(as.character(1:(ncol(A) + 1)))

  if (return.sums) {
    return(colSums(psi))
  } else {
    return(psi)
  }
}


#' MCCS GLM estimating function
#'
#' @param Y outcome, a numeric vector
#' @param A ...
#'
#' @return individual or summation estimating function values
#'
#' @export
get.psi.glm.mccs <- function(Y, Astar, L, g,
                             inv.link, d.inv.link,
                             var.e, B = 10, seed = 123,
                             return.sums = T) {

  # set seed upon each evaluation
  set.seed(seed)

  n <- length(Y)
  len.a <- ncol(Astar)

  # mean of real components of psi0 with simulated imaginary measurement error
  psi <- vapply(
    X = 1:B,
    FUN = function(b) {

      # sample n imaginary measurement error values
      e.tilde <- matrix(complex(
        len.a, real = 0,
        imaginary = mvrnorm(n = n,
                            mu = rep(0, len.a),
                            Sigma = diag(var.e))),
        nrow = n, byrow = F)

      # predictor matrix
      Atilde <- Astar + e.tilde
      Xtilde <- cbind(1, Atilde, L, L * Atilde)

      # evaluate glm psi at Astar + e.tilde and take real component
      Re(get.psi.glm(Y = Y, X = Xtilde, g = g,
                     inv.link = inv.link, d.inv.link = d.inv.link,
                     return.sums = F))
    },

    FUN.VALUE = numeric(n * 2 * (len.a + 1))) |>
    rowMeans() |>
    matrix(nrow = n, ncol = 2 * (len.a + 1), byrow = F)

  if (return.sums) {
    return(colSums(psi))
  } else {
    return(psi)
  }
}


#' MCCS IPW estimating function
#'
#' @param Y outcome, a numeric vector
#' @param A ...
#'
#' @return individual or summation estimating function values
#'
#' @export
get.psi.ipw.mccs <- function(Y, Astar, L, g,
                             inv.link, d.inv.link,
                             var.e, B = 10, seed = 123,
                             return.sums = T, coef.a.l, var.a.l,
                             mean.a, cov.a) {

  # set seed upon each evaluation
  set.seed(seed)

  n <- length(Y)
  len.a <- ncol(Astar)

  # mean of real components of psi0 with simulated imaginary measurement error
  psiMCCS <- vapply(
    X = 1:B,
    FUN = function(b) {

      # sample n imaginary measurement error values
      e.tilde <- matrix(complex(
        len.a, real = 0,
        imaginary = mvrnorm(n = n,
                            mu = rep(0, len.a),
                            Sigma = diag(var.e))),
        nrow = n, byrow = F)

      # evaluate psi0 at Astar + e.tilde and take real component
      Re(get.psi.ipw(Y = Y, L = L, g = g, A = Astar + e.tilde,
                     inv.link = inv.link, d.inv.link = d.inv.link,
                     coef.a.l = coef.a.l, var.a.l = var.a.l,
                     mean.a = mean.a, cov.a = cov.a,
                     return.sums = F)) },

    FUN.VALUE = numeric(n * (len.a + 1))) |>
    rowMeans() |>
    matrix(nrow = n, ncol = len.a + 1, byrow = F)

  if (return.sums) {
    return(colSums((psiMCCS)))
  } else {
    return(psiMCCS)
  }
}






