#' Inverse logit link
#'
#' @param x a numeric vector
#'
#' @return vector of inverse logit values
#'
#' @export
inv.logit <- function(x) {
  exp(x) / (1 + exp(x))
}

#' Derivative of inverse logit
#'
#' @param x
#'
#' @return vector of derivative values
#'
#' @export
d.inv.logit <- function(x) {
  exp(x) / (1 + exp(x)) ^ 2
}

#' Inverse identity link
#'
#' @param x
#'
#' @return x
#'
#' @export
inv.ident <- function(x) {
  x
}

#' Derivative of inverse logit
#'
#' @param x
#'
#' @return a vector of ones
#'
#' @export
d.inv.ident <- function(x) {
  rep(1, length(x))
}

