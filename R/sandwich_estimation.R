#' Sandwich variance estimator
#'
#' @param ghat a numeric vector, estimated parameters
#' @param get.psi estimating function
#' @param n a positive integer, the sample size
#'
#' @return estimated covariance matrix of ghat
#'
#' @export
get.sand.est <- function(ghat, get.psi, n) {

  # D: empirical mean of derivative of Psi
  D <- numDeriv::jacobian(
    f = function(x) colSums(get.psi(x)),
    x = ghat,
    method = "simple") / -n
  Dinv <- solve(D)

  # B: empirical mean of outer product of Psi
  Psi <- get.psi(ghat)
  Omega <- matrix(rowMeans(apply(Psi, 1, function(psi) psi %*% t(psi))),
                  nrow = length(ghat))

  # sandwich estimator
  #sqrt(diag(Dinv %*% Omega %*% t(Dinv) / n)[1:4])
  return(Dinv %*% Omega %*% t(Dinv) / n)
}

#' Bias-corrected sandwich variance estimator
#'
#' @param ghat a numeric vector, estimated parameters
#' @param get.psi estimating function
#' @param n a positive integer, the sample size
#' @param b bias correction parameter
#'
#' @return estimated covariance matrix of ghat
#'
#' @export
get.sand.est.bc <- function(ghat, get.psi, n, b = 0.7) {

  # DD: derivatives of Psi
  DD <- array(numDeriv::jacobian(
    f = function(x) get.psi(x),
    x = ghat,
    method = "simple"),
    dim = c(n, length(ghat), length(ghat)))
  D <- apply(DD, 2:3, sum)
  Dinv <- solve(D)

  # BB: outer products of Psi
  Psi <- get.psi(ghat)
  BB <- simplify2array(lapply(
    1:n, function(i) {
    Psi[i,] %*% t(Psi[i,]) }))

  # HH: diagonals of H matrices for bias correction
  HH <- apply(DD, 1, function(Di){
    (1 - pmin(b, diag(Di %*% Dinv) ) )^(-0.5)
  })

  # BB.bc: bias-corrected B matrices
  BB.bc <- vapply(
    1:n, FUN.VALUE = D,
    FUN = function(i) {
      diag(HH[,i]) %*% BB[,,i] %*% diag(HH[,i])
    })
  B.bc <- apply(simplify2array(BB.bc), 1:2, sum)

  # sandwich estimator
  return(Dinv %*% B.bc %*% t(Dinv))
}


