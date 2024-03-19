#' Create model matrix with complex numbers
#'
#' @inheritParams get.psi.glm
#'
#' @param trms a terms object
#'
#' @return model matrix
#'
#' @export

mod.mat <- function(trms, data) {

  respname <- as.character(attr(trms, "variables")[[attr(trms,
                                                         "response") + 1]])
  termlabs <- attr(trms, "term.labels")
  interactions <- grep(":", attr(trms, "term.labels"), value = TRUE)
  if (length(interactions) == 0) {
    prednames <- termlabs
  } else {
    prednames <- termlabs[!(termlabs %in% interactions)]
  }
  modelframe <- data[prednames]
  if (length(interactions) != 0) {
    for (inter in interactions) {
      intersplit <- strsplit(inter, ":")[[1]]
      modelframe[, inter] <- data[, intersplit[1]] *
        data[, intersplit[2]]
    }
  }
  if (attr(trms, "intercept") == 1) {
    modelmatrix <- as.matrix(data.frame(`(intercept)` = rep(
      1, length(modelframe[, 1])), modelframe))
  } else {
    modelmatrix <- as.matrix(modelframe)
  }
  if (attr(trms, "intercept") == 1) {
    attr(modelmatrix, "assign") <- 0:length(termlabs)
  } else {
    attr(modelmatrix, "assign") <- 1:length(termlabs)
  }
  if (length(interactions) != 0) {
    attr(modelmatrix, "dimnames") <- list(as.character(1:
      length(modelframe[, 1])), c("(intercept)", prednames, interactions))
  } else {
    attr(modelmatrix, "dimnames") <- list(as.character(1:
      length(modelframe[, 1])), c("(intercept)", prednames))
  }
  return(modelmatrix)
}
