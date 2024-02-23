#' Fit oracle doubly robust estimating equation
#'
#' @inheritParams get.psi.glm
#' @inheritParams get.psi.ipw
#'
#' @param start an optional numeric vector, starting parameter values
#' @param a exposure value at which to estimate mean potential outcome
#'
#' @return root of GLM estimating function
#'
#' @export
fit.dr <- function(Y, A, L, a, formula, inv.link, d.inv.link,
                   start = NULL, return.var = TRUE,
                   mean.a = NULL, cov.a = NULL,
                   coef.a.l = NULL, var.a.l = NULL) {

  len.a <- ifelse(is.vector(A), 1, ncol(A))     # dimension of A
  len.msm <- len.a + 1        # dimension of MSM parameters
  len.ps <- 3 * len.a         # dimension of propensity score model parameters
  n <- ifelse(is.vector(A), length(A), nrow(A)) # sample size

  # compute marginal mean and covariance of A if not supplied
  if (is.null(mean.a)) {mean.a <- colMeans(as.matrix(A))}
  if (is.null(cov.a)) {cov.a <- cov(matrix(A))}

  # set starting value if not supplied
  if (is.null(start)) {start <- rep(0, len.msm)}

  # fit propensity score model if not supplied
  if (is.null(coef.a.l) | is.null(var.a.l)) {
    model.a.l <- lm(A ~ L)
    coef.a.l <- t(coef(model.a.l))
    var.a.l <- apply(matrix(model.a.l$residuals), 2, var)
  }

  # solve IPW equation
  root <- tryCatch(
    expr = rootSolve::multiroot(
      f = function(x) {
        get.psi.ipw(
          Y = Y, A = A, L = L, g = x[1:(len.a + 1)],
          inv.link = inv.link, d.inv.link = d.inv.link,
          coef.a.l = coef.a.l, var.a.l = var.a.l,
          mean.a = mean.a, cov.a = cov.a) },
      start = start)$root,
    warning = function(w) {message(w); rep(NA, len.msm)},
    error = function(e) {message(e); rep(NA, len.msm)})

  # estimate E{Y(a)} for each supplied a value
  EYa <- vapply(X = a,
                FUN.VALUE = 0,
                FUN = function(aa) {
                  a.mod.mat <- mod.mat(terms(as.formula(gsub("A", "a", formula))),
                                       data = data.frame(a = aa, L))
                  mean(inv.link(a.mod.mat %*% root))
                })
  names(EYa) <- paste0("EYa.", 1:length(a))
  ghat = c(root, EYa)

  # combine MSM and PS model parameters
  est <- c(root, coef.a.l, log(var.a.l))
  names(est) <- c(
    paste0("g.", 0:len.a),
    paste0("coef.a.l.", 1:(2*len.a)),
    paste0("log.var.a.l", 1:len.a))

  # sandwich variance estimate if requested
  evar = matrix(NA, len.est, len.est)
  if (return.var) {
    evar <- tryCatch(
      expr = get.sand.est(
        ghat = ghat,
        n = n,
        get.psi = function(x) {
          ght <- head(x, -length(a))
          EYa <- tail(x, length(a))
          cbind(
            get.psi.glm(
              Y = Y, A = A, L = L, g = ght, formula = formula,
              inv.link = inv.link, d.inv.link = d.inv.link, return.sums = F),
            vapply(X = 1:length(a),
                   FUN.VALUE = numeric(n),
                   FUN = function(aa) {
                     a.mod.mat <- mod.mat(terms(as.formula(gsub("A", "a", formula))),
                                          data = data.frame(a = a[aa], L))
                     EYa[aa] - inv.link(a.mod.mat %*% ght)
                   }))}),
      warning = function(w) {message(w); matrix(NA, len.est, len.est)},
      error = function(e) {message(e); matrix(NA, len.est, len.est)})
  }

  return(list(est = ghat,
              var = evar))
}
