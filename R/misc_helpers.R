#' Assess estimating equation
#'
#' @param A ...
#'
#' @return individual or summation estimating function values
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
