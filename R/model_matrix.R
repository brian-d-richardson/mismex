#' Create model matrix with complex numbers
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
  interactions <- grep(":", attr(trms, "term.labels"), value = TRUE)

  # Handle polynomial terms of the form "I(A^n)"
  polynomial_terms <- grep("^I\\(A\\^.*\\)$", termlabs, value = TRUE)

  # Non-interaction, non-polynomial terms
  non_polynomial_terms <- termlabs[!(termlabs %in% c(interactions, polynomial_terms))]

  # Check for "L" and separate it to handle ordering
  L_terms <- grep("^L$", non_polynomial_terms, value = TRUE)
  non_polynomial_terms <- setdiff(non_polynomial_terms, L_terms)

  # Initialize model frame with non-polynomial, non-L terms
  modelframe <- data[non_polynomial_terms]

  # Add polynomial terms for A before L
  if (length(polynomial_terms) != 0) {
    for (poly_term in polynomial_terms) {
      # Extract the variable and exponent from terms like "I(A^2)"
      term_split <- strsplit(poly_term, "\\^|\\(|\\)")[[1]]
      var_name <- term_split[2]  # The variable name (A)
      exponent <- as.numeric(term_split[3])  # The exponent
      # Compute the polynomial term
      modelframe[, poly_term] <- data[, var_name]^exponent
    }
  }

  # Add L terms after polynomial terms
  if (length(L_terms) != 0) {
    modelframe[, L_terms] <- data[, L_terms]
  }

  # Add interaction terms to model frame
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

  # Set column names with the correct ordering: (Intercept), non-polynomial terms, polynomials, L, interactions
  colnames <- c("(Intercept)", non_polynomial_terms, polynomial_terms, L_terms, interactions)
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
