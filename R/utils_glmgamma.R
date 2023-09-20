#' Compute score function for a generalized linear model with Gamma response.
#'
#' @param fit TBD
#'
#' @param theta a numeric vector of length (p+1), containing MLE of parameters in a linear model.
#'
#' @return Score matrix with n rows and (p+2) columns.
#'
glmScorePIT = function(fit, theta){

  # MLE of shape parameter
  mle.alpha<- theta

  # Extract family function, design matrix (x) from fit object, compute the linear predictor and extract response variable.
  fm              <- family(fit)
  x               <- fit$x
  y               <- fit$y
  mle.coef        <- coef(fit)
  linearPredictor <- x %*% mle.coef

  # Check if we need to add offset to the linear predictor.
  if( !is.null(fit$offset) ){
    linearPredictor <- linearPredictor + fit$offset
  }

  # Calculate muhat, estimated value of mean function in GLM by calling linkinv
  muhat          <- fm$linkinv(linearPredictor)

  # Compute the score matrix.
  # Family function from stats package to get the link function and the
  # derivative of the inverse-link function with respect to the linear predictor
  S1    <- 1 - digamma(mle.alpha) + log(mle.alpha) + log(y/muhat) - (y/muhat)
  S2    <- ( fm$mu.eta(linearPredictor) / muhat ) * mle.alpha * ( (y/muhat) - 1 )
  S2    <- as.vector(S2) * x
  S     <- cbind(S1,S2)

  # Calculate probability inverse transfer of data
  pit <- pgamma(q = as.numeric(y), shape = mle.alpha, scale = muhat/mle.alpha)

  # Return the score function and pit values
  return( list(Score = S, pit = pit) )

}


#' Compute maximum likelihood estimates for a generalized linear model with Gamma response.
#'
#' @param fit is an object of class \code{glm} and its default value is NULL. If a fit of class \code{glm} is provided,
#' the arguments \code{x}, \code{y}, and \code{l} will be ignored. We recommend using \code{\link{glm2}} function from
#' \code{\link{glm2}} package since it provides better convergence while optimizing the likelihood to estimate
#' coefficients of the model by IWLS method. It is required to return design matrix by \code{x} = \code{TRUE} in
#' \code{\link{glm}} or \code{\link{glm2}} function. For more information on how to do this, refer to the help
#' documentation for the \code{\link{glm}} or \code{\link{glm2}} function.
#'
#' @return a numeric vector of estimates.
#'
glmMLE = function(fit){

  # coef function from stats package extracts MLE estimate of beta parameter in the model.
  mle.coef                 <- coef(fit)

  # Use gamma.shape function from MASS package to compute MLE of shape parameter in Gamma distribution
  mle.alpha                <- MASS::gamma.shape(fit, it.lim = 30, eps.max = 1e-9)[1]$alpha

  # Build a vector for MLE estimates, p parameters for beta and last parameter for mle of shape parameter
  par                      <- c(mle.coef, mle.alpha)
  names(par)[length(par)]  <- 'shape'

  return(par)

}


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
glmgammaFisherByHessian = function(fit, mle_shape){

  # Extract MLE of coefficients of hte model from fit object
  mle.coef <- coef(fit)

  # Extract the design matrix
  X        <- fit$x
  y        <- fit$y
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
  Xmm             <- as.vector( sqrt(y/miohat) ) * X


  # Compute Hessian matrix
  p <- ncol(X) + 1
  Hessian     <- matrix(0, nrow = p, ncol = p)

  first_term  <- -mle_shape * t(fit$x) %*% fit$x / n

  second_term <- mle_shape * t(Xm) %*% Xm / n

  third_term  <- mle_shape * t(Xmm) %*% ( Xmm ) / n

  fourth_term <- -2*mle_shape * t(Xmm) %*% Xmm / n

  Hessian[1,1]     <- (-1)*trigamma(mle_shape) + (1)/(mle_shape)
  Hessian[2:p,2:p] <- first_term + second_term + third_term + fourth_term

  return(Hessian)

}
