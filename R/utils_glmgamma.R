#' Compute Fisher information matrix in the case of GLM with Gamma response variable.
#'
#' @description Computes the Fisher information matrix of parameters in the case of a generalized linear model with the assumption of Gamma
#' for the response variable.
#'
#' @param par mle of shape parameter for the Gamma distribution.
#'
#' @param x matrix of explanatory variable.
#'
#' @return a Fisher information matrix with n rows (number of observations in the sample) and p columns (number of parameters in the model).
#'
#' @noRd
observedHessianMatrixGLMGamma = function(par, x){

  n           <- nrow(x)
  rindx       <- ncol(x) + 1
  FI          <- matrix(NA, nrow = ncol(x)+2, ncol = ncol(x)+2)
  FI[1,1]     <- trigamma(par) - (1/par)
  FI[2:rindx,2:rindx] <- (-par/n) * ( t(X) %*% X )
  FI[1,2:rindx]   <- FI[2:rindx,1] <- 0

  return(FI)

}
