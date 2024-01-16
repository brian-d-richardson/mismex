#' Fit oracle GLM estimating equation
#'
#' @param Y outcome, a numeric vector
#' @param A ...
#'
#' @return individual or summation estimating function values
#'
#' @export
fit.glm <- function(Y, X, inv.link, d.inv.link, g.start = NULL) {

  if (is.null(g.start)) {
    g.start <- rep(0, ncol(X) + 1)
  }

  # Solve oracle IPW equation
  root <- rootSolve::multiroot(
    f = function(g)  get.psi.glm(g = g, Y = Y, X = X,
                                 inv.link = inv.link, d.inv.link = d.inv.link),
    start = rep(0, ncol(X)))

  return(root$root)
}
#fit.glm(Y = Y, X = cbind(1, A, L, A * L),
#        inv.link = inv.logit, d.inv.link = d.inv.logit)


#' Fit oracle IPW estimating equation
#'
#' @param Y outcome, a numeric vector
#' @param A ...
#'
#' @return individual or summation estimating function values
#'
#' @export
fit.ipw <- function(Y, A, L,
                    inv.link, d.inv.link,
                    g.start = NULL) {

  len.a <- ncol(A)
  mean.a <- colMeans(A)
  cov.a <- cov(A)
  if (is.null(g.start)) {
    g.start <- rep(0, len.a + 1)
  }

  # get starting values for ps model
  model.a.l <- lm(A ~ L)
  start <- c(g.start,
             as.numeric(t(coef(model.a.l))),
             log(apply(model.a.l$residuals, 2, var)))

  # Solve oracle IPW equation
  root <- rootSolve::multiroot(
    f = function(x) {

        psi.ps <- get.psi.ps(
          A = A, L = L,
          coef.a.l = matrix(x[len.a + 1 + 1:(2 * len.a)], nrow = len.a),
          var.a.l = exp(x[3 * len.a + 1 + 1:len.a]))

        psi.ipw <- get.psi.ipw(
          Y = Y, A = A, L = L, g = x[1:(len.a + 1)],
          inv.link = inv.link, d.inv.link = d.inv.link,
          coef.a.l = matrix(x[len.a + 1 + 1:(2 * len.a)], nrow = len.a),
          var.a.l = exp(x[3 * len.a + 1 + 1:len.a]),
          mean.a = mean.a, cov.a = cov.a)

      return(c(psi.ps, psi.ipw)) },
    start = start)

  ret <- root$root
  names(ret) <- c(
    paste0("g.", 0:len.a),
    paste0("coef.a.l.", 1:(2 * len.a)),
    paste0("log.var.a.l.", 1:len.a))

  return(ret)
}
#fit.ipw(Y = Y, A = A, L = L,
#        inv.link = inv.logit, d.inv.link = d.inv.logit,
#        g.start = coef(glm(Y ~ A, family = binomial)))


#' Fit MCCS GLM estimating equation
#'
#' @param Y outcome, a numeric vector
#' @param A ...
#'
#' @return individual or summation estimating function values
#'
#' @export
fit.glm.mccs <- function(Y, Astar, L,
                         var.e, B = 10, seed = 123,
                         inv.link, d.inv.link,
                         g.start = NULL) {

  len.a <- ncol(Astar)

  # get naive estimates to use as starting values
  root.naive <- fit.glm(Y = Y, X = cbind(1, Astar, L, Astar * L),
                        inv.link = inv.link, d.inv.link = d.inv.link,
                        g.start = g.start)

  # Solve MCCS GLM equation
  root <- rootSolve::multiroot(
    f = function(x) {

      # IPW equation
      get.psi.glm.mccs(
        Y = Y, Astar = Astar, L = L, g = x,
        inv.link = inv.link, d.inv.link = d.inv.link,
        var.e = var.e, B = B, seed = seed) },

    start = root.naive)

  ret <- root$root
  names(ret) <- paste0("g.", 0:(length(ret) - 1))
  return(ret)
}
#fit.glm.mccs(Y = Y, Astar = Astar, L = L,
#             var.e = var.e, B = 10, seed = 123,
#             inv.link = inv.logit, d.inv.link = d.inv.logit)


#' Fit MCCS IPW estimating equation
#'
#' @param Y outcome, a numeric vector
#' @param A ...
#'
#' @return individual or summation estimating function values
#'
#' @export
fit.ipw.mccs <- function(Y, Astar, L,
                         var.e, B = 10, seed = 123,
                         inv.link, d.inv.link,
                         g.start = NULL) {

  len.a <- ncol(Astar)
  mean.a <- colMeans(Astar)
  cov.a <- cov(Astar) - diag(var.e)

  # get naive estimates to use as starting values
  root.naive <- fit.ipw(Y = Y, A = Astar, L = L,
                        inv.link = inv.link, d.inv.link = d.inv.link,
                        g.start = g.start)

  # Solve MCCS IPW equation
  root <- rootSolve::multiroot(
    f = function(x) {

      # probability weight equation
      psi.ps <- get.psi.ps(
        A = Astar, L = L,
        coef.a.l = matrix(x[len.a + 1 + 1:(2 * len.a)], nrow = len.a),
        var.a.l = exp(x[3 * len.a + 1 + 1:len.a]))

      # IPW equation
      psi.ipw <- get.psi.ipw.mccs(
        Y = Y, Astar = Astar, L = L, g = x[1:(len.a + 1)],
        var.e = var.e, B = B, seed = seed,
        inv.link = inv.link, d.inv.link = d.inv.link,
        coef.a.l = matrix(x[len.a + 1 + 1:(2 * len.a)], nrow = len.a),
        var.a.l = exp(x[3 * len.a + 1 + 1:len.a]),
        mean.a = mean.a, cov.a = cov.a)

      return(c(psi.ps, psi.ipw)) },
    start = root.naive)

  ret <- root$root
  names(ret) <- c(
    paste0("g.", 0:len.a),
    paste0("coef.a.l.", 1:(2 * len.a)),
    paste0("log.var.a.l.", 1:len.a))
  return(ret)
}
#fit.ipw.mccs(Y = Y, Astar = Astar, L = L,
#             var.e = var.e, B = 10, seed = 123,
#             inv.link = inv.logit, d.inv.link = d.inv.logit)
















