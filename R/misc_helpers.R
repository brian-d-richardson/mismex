#' Assess estimating equation
#'
#' @param ee a matrix of numbers, columns correspond to components of the estimating equation, rows correspond to iid replicates
#' @param digits a non-negative integer, the number of digits displayed in return values (default is 2)
#'
#'
#'#' @return a list with the following elements:
#' \itemize{
#' \item{means: mean values of the estimating function}
#' \item{t.stats: t-statistics to test the hypotheses of mean zero components}
#' }
#' @return a named list with mean values and t-statistics to test the hypothesis of mean zero
#'
#' @export
assess.ee <- function(ee, digits = 2) {

  # means of estimating equation values (for each component)
  means <- colMeans(ee)

  # t statistics of H_0: mean = 0 (for each component)
  n <- nrow(ee)
  t.stats <- apply(ee, 2, function(x) mean(x) * sqrt(n) / sd(x))

  return(list(means = round(means, digits),
              t.stats = round(t.stats, digits)))
}

# Format data for dose response curve plot
#'
#' @param a a matrix of numbers, columns correspond to components of the estimating equation, rows correspond to iid replicates
#'
#' @return a data frame
#'
#' @export
format.gfmla.res <- function(a, res, alpha = 0.05) {
  EYa <- tail(res$est, length(a))
  se <- sqrt(tail(diag(res$var), length(a)))
  drc.dat <- data.frame(
    a = a,
    est = EYa,
    se = se,
    lower = EYa - qnorm(1 - alpha / 2) * se,
    upper = EYa + qnorm(1 - alpha / 2) * se)
  return(drc.dat)
}
