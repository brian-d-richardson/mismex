#' Sandwich variance estimator
#'
#' @inheritParams get.psi.ipw
#'
#' @param ghat a numeric vector, estimated parameters
#' @param get.psi estimating function
#'
#' @return estimated covariance matrix of ghat
#'
#' @export
get.sand.est <- function(ghat, get.psi, n = NULL) {

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

