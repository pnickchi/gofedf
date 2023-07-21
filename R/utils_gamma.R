#' Compute maximum likelihood estimates of shape and scale parameter in Gamma distribution.
#'
#' @description Estimate the shape and scale parameters of the Gamma distribution using the method of maximum likelihood,
#' and apply the Newton-Raphson method on the profile log-likelihood to estimate the shape parameter.
#'
#' @param x a numeric vector.
#'
#' @param ur logical. If \code{TRUE} the rate parameter is returned. Otherwise the scale is returned.
#'
#' @return a vector of length two with shape and scale/rate.
#'
#' @noRd
estimateGamma <- function(x, ur){

  # Find the the number of observations in the sample
  n <- length(x)

  # Compute mean and variance of the sample
  m <- mean(x)
  s <- var(x)

  # Initialize values
  b <- s / m
  a <- m / b
  mlog <- mean(log(x))
  logm <- log(m)
  shape_old <- a
  shape_new <- shape_old -(log(shape_old)-logm + mlog -digamma(shape_old))/(1/shape_old-trigamma(shape_old))
  bnew <- m/shape_new

  # Check if the new shape is negative and replace accordingly
  if( shape_new < 0) shape_new <- shape_old/2

  # Optimization part
  while ( abs(shape_new-shape_old) > 1e-7){
    shape_old <- shape_new
    old.score = (log(shape_old)-log(m)+ mlog -digamma(shape_old))
    old.score.derivative = 1/shape_old-trigamma(shape_old)
    shape_new <- shape_old - old.score/old.score.derivative
    if( shape_new < 0) shape_new <- shape_old/2
  }

  beta_val  <- m / shape_new
  alpha     <- shape_new

  # Return MLE estimates
  if(ur){
    return(c(alpha,1/beta_val))
  }else{
    return(c(alpha, beta_val))
  }

}
#' Compute observed Hessian matrix evaluated at the sample for Gamma distribution.
#'
#' @param par a numeric vector with maximum likelihood estimates of shape and scale.
#'
#' @return Hessian matrix
#'
#' @noRd
observedHessianMatrixGamma = function(par){

  alpha.hat  <- par[1]
  lambda.hat <- par[2]

  res      <- matrix(0, nrow = 2, ncol = 2)
  res[1,1] <- trigamma(alpha.hat)
  res[1,2] <- res[2,1] <- -1/lambda.hat
  res[2,2] <- alpha.hat / lambda.hat^2
  return(res)

}

# An internal function to compute F(x) function for Gamma distribution. Used for example.
GammaFx = function(x, theta){

  if( x > 0){
    res <- pgamma(q = x, shape = theta[1], rate = theta[2])
  }else{
    res <- 0
  }

  return(res)
}

# An internal function to compute score function for Gamma distribution. Used for example.
GammaScore = function(x, theta){

  alpha  <- theta[1]
  lambda <- theta[2]

  S1 <- log(lambda) - digamma(alpha) + log(x)
  S2 <- (alpha / lambda) - x
  S  <- cbind(S1,S2)

  return(S)
}
