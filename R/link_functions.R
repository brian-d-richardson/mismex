#' Inverse logit
#'
#' @param x
#'
#' @return
#'
#' @export
inv.logit <- function(x) {
  exp(x) / (1 + exp(x))
}

#' Derivative of inverse logit
#'
#' @param x
#'
#' @return
#'
#' @export
d.inv.logit <- function(x) {
  exp(x) / (1 + exp(x)) ^ 2
}

#' indentity
#'
#' @param x
#'
#' @return
#'
#' @export
inv.ident <- function(x) {
  x
}

#' Derivative of inverse logit
#'
#' @param x
#'
#' @return
#'
#' @export
d.inv.ident <- function(x) {
  1
}

