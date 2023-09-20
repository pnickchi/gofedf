#' Compute score function for Gamma distribution.
#'
#' @param x a numeric vector of length n
#'
#' @param theta a numeric vector of length two, containing MLE of parameters in Gamma dist.
#'
#' @return Score matrix with n rows and two columns.
#'
gammaScore = function(x, theta){

  # Extract MLE of shape and scale from theta argument
  alpha  <- theta[1]
  lambda <- theta[2]

  # Compute the first and second column of score matrix
  S1 <- log(lambda) - digamma(alpha) + log(x)
  S2 <- (alpha / lambda) - x

  # Create the score matrix with n rows and two columns
  S  <- cbind(S1,S2)

  # Return score matrix
  return(S)
}


#' Compute maximum likelihood estimate of shape and scale parameter in Gamma distribution.
#'
#' @description Estimate the MLE of shape and scale parameters of the Gamma distribution using the
#' Newton-Raphson method on the profile log-likelihood to estimate the shape parameter.
#'
#' @param x a numeric vector of length n
#'
#' @param ur logical. If \code{TRUE} the rate parameter is returned. Otherwise the scale is returned.
#'
#' @return a vector of length two with shape and scale/rate.
#'
gammaMLE = function(x, ur){

  # Find the the number of observations in the sample
  n <- length(x)

  # Compute mean and variance of the sample
  m <- mean(x)
  s <- var(x)

  # Initialize values
  b         <- s / m
  a         <- m / b
  mlog      <- mean(log(x))
  logm      <- log(m)
  shape_old <- a
  shape_new <- shape_old -(log(shape_old)-logm + mlog -digamma(shape_old))/(1/shape_old-trigamma(shape_old))
  bnew      <- m/shape_new

  # Check if the new shape is negative and replace accordingly
  if( shape_new < 0) shape_new <- shape_old/2

  # Optimization part
  while ( abs(shape_new-shape_old) > 1e-7){
    shape_old <- shape_new
    old.score <- (log(shape_old)-log(m)+ mlog -digamma(shape_old))
    old.score.derivative <- 1/shape_old-trigamma(shape_old)
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


#' Compute probability inverse transform values for Gamma distribution
#'
#' @param x a numeric vector of length n
#'
#' @param theta a numeric vector of length two, containing MLE of parameters in Gamma dist.
#'
#' @return a vector of length n containing the probability inverse transformed (PIT) values
#'
gammaPIT = function(x, theta){

  # Extract MLE of shape and scale from theta argument
  a <- theta[1]
  l <- theta[2]

  # Compute PIT values
  res <- pgamma(q = x, shape = a, rate = l)

  return(res)
}


#' Compute Fisher information matrix by the negative expected value of Hessian matrix in Gamma distribution.
#'
#' @param theta a numeric vector of length two, containing MLE of parameters in Gamma dist
#'
#' @return Fisher information matrix for Gamma distribution
#'
gammaFisherByHessian = function(theta){

  # Extract MLE of shape and scale from theta argument
  alpha.hat  <- theta[1]
  lambda.hat <- theta[2]

  # Define a 2x2 matrix
  res      <- matrix(0, nrow = 2, ncol = 2)

  # Compute values
  res[1,1] <- trigamma(alpha.hat)
  res[1,2] <- res[2,1] <- -1/lambda.hat
  res[2,2] <- alpha.hat / lambda.hat^2

  # Return the Fisher information matrix
  return(res)

}

