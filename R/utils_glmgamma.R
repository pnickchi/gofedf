#' Compute Fisher information matrix in the case of GLM with Gamma response variable.
#'
#' @description Computes the Fisher information matrix of parameters in the case of a generalized linear model with the assumption of Gamma
#' for the response variable.
#'
#' @param mle_shape mle of shape parameter for the Gamma distribution.
#'
#' @param fit is an object of class \code{glm} and its default value is NULL. If a fit of class \code{glm} is provided,
#' the arguments \code{x}, \code{y}, and \code{l} will be ignored. We recommend using \code{\link{glm2}} function from
#' \code{\link{glm2}} package since it provides better convergence while optimizing the likelihood to estimate
#' coefficients of the model by IWLS method. It is required to return design matrix by \code{x} = \code{TRUE} in
#' \code{\link{glm}} or \code{\link{glm2}} function. For more information on how to do this, refer to the help
#' documentation for the \code{\link{glm}} or \code{\link{glm2}} function.
#'
#' @return a Fisher information matrix with n rows (number of observations in the sample) and p columns (number of parameters in the model).
#'
#' @noRd
observedHessianMatrixGLMGamma = function(fit, mle_shape){

  # Extract MLE of coefficients of hte model from fit object
  mle.coef <- coef(fit)

  # Extract the design matrix
  X        <- fit$x
  n        <- nrow(X)

  # Compute the linear predictor value
  linearPredictor <- X %*% mle.coef

  # Check if we need to add any offset to the linear predictor
  if( !is.null(fit$offset) ){
    linearPredictor <- linearPredictor + fit$offset
  }

  # Extarct the family function from fit object
  fm       <- family(fit)

  # Calculate miohat, estimated value of mean function in GLM
  miohat          <- fm$linkinv(linearPredictor)
  partial.mu      <- fm$mu.eta(linearPredictor)
  temp            <- partial.mu / miohat
  Xm              <- X * as.vector(temp)

  # Compute Hessian matrix
  p <- ncol(X) + 1
  Hessian          <- matrix(0, nrow = p, ncol = p)
  Hessian[1,1]     <- trigamma(mle_shape) - (1/mle_shape)
  Hessian[2:p,2:p] <- mle_shape * t(Xm) %*% Xm / n

  return(Hessian)

}
