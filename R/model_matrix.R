#' Create model matrix for A|L model with complex numbers
#'
#' @inheritParams get.psi.glm
#'
#' @param trms a terms object
#'
#' @return model matrix
#'
#' @export
#'
#'
mod.mat <- function(trms, data) {

  respname <- as.character(attr(trms, "variables")[[attr(trms, "response") + 1]])
  termlabs <- attr(trms, "term.labels")
  interactions <- grep(":", termlabs, value = TRUE)

  # Handle polynomial terms for A* and L* (A, A1, A2, ..., and L, L1, L2, ...)
  A_polynomial_terms <- grep("^I\\(A[0-9]*\\^.*\\)$", termlabs, value = TRUE)
  L_polynomial_terms <- grep("^I\\(L[0-9]*\\^.*\\)$", termlabs, value = TRUE)

  # Identify non-polynomial, non-interaction terms
  non_polynomial_terms <- setdiff(termlabs, c(interactions, A_polynomial_terms, L_polynomial_terms))

  # Separate plain "A*" and "L*" terms (e.g., A, A1, L, L1)
  A_terms <- grep("^A[0-9]*$", non_polynomial_terms, value = TRUE)
  L_terms <- grep("^L[0-9]*$", non_polynomial_terms, value = TRUE)
  non_polynomial_terms <- setdiff(non_polynomial_terms, c(A_terms, L_terms))

  # Initialize model frame with non-polynomial, non-A, non-L terms
  modelframe <- data[non_polynomial_terms]

  # Add "A*" terms and A polynomials
  if (length(A_terms) != 0) {
    for (term in A_terms) {
      modelframe[, term] <- data[, term]
    }
  }
  if (length(A_polynomial_terms) != 0) {
    for (term in A_polynomial_terms) {
      term_split <- strsplit(term, "\\^|\\(|\\)")[[1]]
      var_name <- term_split[2]  # Variable name (e.g., A, A1)
      exponent <- as.numeric(term_split[3])  # Extract exponent for A
      modelframe[, term] <- data[, var_name]^exponent
    }
  }

  # Add "L*" terms and L polynomials
  if (length(L_terms) != 0) {
    for (term in L_terms) {
      modelframe[, term] <- data[, term]
    }
  }
  if (length(L_polynomial_terms) != 0) {
    for (term in L_polynomial_terms) {
      term_split <- strsplit(term, "\\^|\\(|\\)")[[1]]
      var_name <- term_split[2]  # Variable name (e.g., L, L1)
      exponent <- as.numeric(term_split[3])  # Extract exponent for L
      modelframe[, term] <- data[, var_name]^exponent
    }
  }

  # Add interaction terms
  if (length(interactions) != 0) {
    for (inter in interactions) {
      intersplit <- strsplit(inter, ":")[[1]]
      modelframe[, inter] <- data[, intersplit[1]] * data[, intersplit[2]]
    }
  }

  # Handle the intercept
  if (attr(trms, "intercept") == 1) {
    modelmatrix <- as.matrix(data.frame(`(Intercept)` = rep(1, nrow(modelframe)), modelframe))
  } else {
    modelmatrix <- as.matrix(modelframe)
  }

  # Assign attributes to the model matrix
  if (attr(trms, "intercept") == 1) {
    attr(modelmatrix, "assign") <- 0:length(termlabs)
  } else {
    attr(modelmatrix, "assign") <- 1:length(termlabs)
  }

  # Set column names with specified order: Intercept, A terms, A polynomials, L terms, L polynomials, interactions
  colnames <- c("(Intercept)", A_terms, A_polynomial_terms, L_terms, L_polynomial_terms, interactions)
  attr(modelmatrix, "dimnames") <- list(as.character(1:nrow(modelframe)), colnames)

  return(modelmatrix)
}



mod.mat.old <- function(trms, data) {

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
