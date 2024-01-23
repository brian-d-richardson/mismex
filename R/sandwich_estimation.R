#' Sandwich variance estimator
#'
#' @inheritParams get.psi.ipw
#'
#' @param ghat a numeric vector, estimated parameters
#' @param oracle logical, indicator for whether true exposure A is provided;
#'    when oracle = TRUE var.e, B, and seed need not be supplied
#'
#' @return estimated covariance matrix of ghat
#'
#' @export
get.sand.est <- function(Y, A, L, ghat, method, oracle = T,
                         inv.link, d.inv.link,
                         var.e = NULL, B = NULL, seed = NULL) {

  len.a <- ncol(A)
  n <- length(Y)

  # if oracle, then set measurement error to 0 and MC replicates to 1
  if (oracle) {
    var.e <- rep(0, len.a)
    B <- 1
    seed <- 1
  }

  # use given method to choose estimating function
  if (method == "GLM") {

    get.psi <- function(g, return.sums) get.psi.glm.mccs(
      Y = Y, Astar = A, L = L, g = g,
      var.e = var.e, B = B, seed = seed,
      inv.link = inv.link, d.inv.link = d.inv.link,
      return.sums = return.sums)

  } else if (method == "IPW") {

    get.psi <- function(g, return.sums) {

      psi <- cbind(

        get.psi.ps(
          A = A, L = L,
          coef.a.l = matrix(g[len.a + 1 + 1:(2 * len.a)], nrow = len.a, byrow = F),
          var.a.l = exp(g[3 * len.a + 1 + 1:len.a]) + var.e,
          return.sums = F),

        get.psi.ipw.mccs(
          Y = Y, Astar = A, L = L, g = g[1:((len.a) + 1)],
          var.e = var.e, B = B, seed = seed,
          inv.link = inv.link, d.inv.link = d.inv.link,
          coef.a.l = matrix(g[len.a + 1 + 1:(2 * len.a)], nrow = len.a, byrow = F),
          var.a.l = exp(g[3 * len.a + 1 + 1:len.a]),
          mean.a = mean(A), cov.a = cov(A) - diag(var.e),
          return.sums = F))

      if (return.sums) {
        return(colSums(psi))
      } else {
        return(psi)
      }
    }
  } else {
    stop("invalid method (expecting GLM or IPW)")
  }

  # D: empirical mean of derivative of Psi
  D <- numDeriv::jacobian(
    f = function(x) get.psi(x, return.sums = T),
    x = ghat,
    method = "simple") / -n
  Dinv <- solve(D)

  # B: empirical mean of outer product of Psi
  Psi <- get.psi(ghat, return.sums = F)
  Omega <- matrix(rowMeans(apply(Psi, 1, function(psi) psi %*% t(psi))),
                  nrow = length(ghat))

  # sandwich estimator
  #sqrt(diag(Dinv %*% Omega %*% t(Dinv) / n)[1:4])
  return(Dinv %*% Omega %*% t(Dinv) / n)
}
