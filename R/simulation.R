#' run one IPW simulation
#'
#' @param A ...
#'
#' @return individual or summation estimating function values
#'
#' @export
sim1 <- function(n = 2000,
                 B = 10,
                 seed = 1) {

  ## for troubleshooting
  #library(devtools); load_all()
  #n = 2000; B = 10; seed = 1;
  
  ## define parameters
  gg = c(-.5, .5, .5, -.5);                         # MSM parameters
  inv.link = function(x) 0.5 * inv.logit(x);        # MSM link function
  d.inv.link = function(x) 0.5 * d.inv.logit(x)     # MSM derivative of link
  var.e <- c(1, 1, 0)                               # measurement error variance
  coef.a.l <- matrix(data = c(1, 1, 0.5, 0, 0, -1), # coefs in A|L model
                     nrow = 3, byrow = T)         
  var.a.l <- c(1, 1, 1)                             # variance of A|L
  
  ## generate data
  set.seed(seed)                                   # seed for reproducibility
  L <- runif(n)                                    # confounder
  A <- mvrnorm(n = n,                              # true exposure                
               mu = c(0, 0, 0),
               Sigma = diag(var.a.l)) +
    cbind(1, L) %*% t(coef.a.l)
  Astar <- A + mvrnorm(n = n,                      # mismeasured exposure
                       m = c(0, 0, 0),
                       Sigma = diag(var.e))
  Y_prob <- L * inv.logit(cbind(1, A) %*% gg)      # mean of binary outcome
  Y <- rbinom(n, 1, Y_prob)                        # binary outcome

  ## estimate MSM parameters
  # (i) oracle logistic regression
  ghat.OL <- tryCatch(
    expr = fit.glm(Y = Y, X = cbind(1, A, L, A*L),
                   inv.link = inv.logit, d.inv.link = d.inv.logit),
    warning = function(w) rep(NA, 8),
    error = function(e) rep(NA, 8))

  #evar.OL <- tryCatch(
  #  expr = get.sand.est(Y = Y, A = A, L = L, ghat = ghat.OL,
  #                      method = "GLM", oracle = T,
  #                      inv.link = inv.logit, d.inv.link = d.inv.logit),
  #  warning = function(w) matrix(NA, 8, 8),
  #  error = function(e) matrix(NA, 8, 8))

  # (ii) naive logistic regression
  ghat.NL <- tryCatch(
    expr = fit.glm(Y = Y, X = cbind(1, Astar, L, Astar*L),
                   inv.link = inv.logit, d.inv.link = d.inv.logit),
    warning = function(w) rep(NA, 8),
    error = function(e) rep(NA, 8))

  #evar.NL <- tryCatch(
  #  expr = get.sand.est(Y = Y, A = Astar, L = L, ghat = ghat.NL,
  #                      method = "GLM", oracle = T,
  #                      inv.link = inv.logit, d.inv.link = d.inv.logit),
  #  warning = function(w) matrix(NA, 8, 8),
  #  error = function(e) matrix(NA, 8, 8))

  # (iii) MCCS logistic regression
  ghat.CL <- tryCatch(
    expr = fit.glm.mccs(Y = Y, Astar = Astar, L = L, var.e = var.e,
                        inv.link = inv.logit, d.inv.link = d.inv.logit,
                        B = B, seed = 123),
    warning = function(w) rep(NA, 8),
    error = function(e) rep(NA, 8))

  #evar.CL <- tryCatch(
  #  expr = get.sand.est(Y = Y, A = A, L = L, ghat = ghat.CL,
  #                      method = "GLM", oracle = T,
  #                      inv.link = inv.logit, d.inv.link = d.inv.logit),
  #  warning = function(w) matrix(NA, 8, 8),
  #  error = function(e) matrix(NA, 8, 8))

  # (iv) oracle IPW estimator
  ghat.OI <- tryCatch(
    expr = fit.ipw(Y = Y, A = A, L = L,
                   inv.link = inv.link, d.inv.link = d.inv.link),
    warning = function(w) rep(NA, 13),
    error = function(e) rep(NA, 13))

  #evar.OI <- tryCatch(
  #  expr = get.sand.est(Y = Y, A = A, L = L, ghat = ghat.OI,
  #                      method = "IPW", oracle = T,
  #                      inv.link = inv.logit, d.inv.link = d.inv.logit),
  #  warning = function(w) matrix(NA, 13, 13),
  #  error = function(e) matrix(NA, 13, 13))

  # (iv) oracle IPW estimator
  ghat.NI <- tryCatch(
    expr = fit.ipw(Y = Y, A = Astar, L = L,
                   inv.link = inv.link, d.inv.link = d.inv.link),
    warning = function(w) rep(NA, 13),
    error = function(e) rep(NA, 13))

  #evar.NI <- tryCatch(
  #  expr = get.sand.est(Y = Y, A = Astar, L = L, ghat = ghat.NI,
  #                      method = "IPW", oracle = T,
  #                      inv.link = inv.logit, d.inv.link = d.inv.logit),
  #  warning = function(w) matrix(NA, 13, 13),
  #  error = function(e) matrix(NA, 13, 13))

  # (vi) MCCS IPW estimator
  ghat.CI <- tryCatch(
    expr = fit.ipw.mccs(Y = Y, Astar = Astar, L = L,
                        var.e = var.e, B = B, seed = 123,
                        inv.link = inv.link, d.inv.link = d.inv.link),
    warning = function(w) rep(NA, 13),
    error = function(e) rep(NA, 13))

  #evar.CI <- tryCatch(
  #  expr = get.sand.est(Y = Y, A = Astar, L = L, ghat = ghat.CI,
  #                      method = "IPW", oracle = F,
  #                      var.e = var.e, B = B, seed = 123,
  #                      inv.link = inv.logit, d.inv.link = d.inv.logit),
  #  warning = function(w) matrix(NA, 13, 13),
  #  error = function(e) matrix(NA, 13, 13))

  ret <- c(n, B, seed,
           ghat.OL[1:4], ghat.NL[1:4], ghat.CL[1:4],
           ghat.OI[1:4], ghat.NI[1:4], ghat.CI[1:4])#,
           #diag(evar.OL)[1:4], diag(evar.NL)[1:4], diag(evar.CL)[1:4],
           #diag(evar.OI)[1:4], diag(evar.NI)[1:4], diag(evar.CI)[1:4])

  names(ret) <- c(
    "n", "B", "seed",
    apply(tidyr::expand_grid(
      "ghat",
      #c("ghat", "evar"),
      c("OL", "NL", "CL", "OI", "NI", "CI"),
      1:4), 1, paste, collapse="."))

  return(ret)
}

#library(tictoc); tic("one sim"); sim.res <- sim1(); toc()
#round(sim.res, 2)
