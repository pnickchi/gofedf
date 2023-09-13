#' Compute score function for Normal dist
#'
#' @param x a numeric vector of length n
#'
#' @param theta a numeric vector of length two, containing MLE of parameters in Normal dist.
#'
#' @return Score matrix with n rows and two columns.
#'
normalScore = function(x, theta){

  # Get the sample size
  n <- length(x)

  # Extract the MLE of mean and sd from theta argument
  m <- theta[1]
  s <- theta[2]

  # Compute the first and second column of score matrix
  S1 <- (x - m) / (s^2)
  S2 <- (x - m)^2/s^3 - rep(1/s, n)

  # Create the score matrix with n rows and two columns
  S  <- cbind(S1,S2)

  # Return score matrix
  return(S)
}


#' Compute MLE estimate for Normal
#'
#' @param x a numeric vector of length n
#'
#' @return a numeric vector of length two, containing MLE of parameters in Normal dist.
#'
normalMLE = function(x){

  # Get the sample size
  n <- length(x)

  # Compute MLE of mean and sd
  m <- mean(x)
  s <- sqrt( (n-1)* var(x) / n)

  # Create a vector of MLE
  theta <- c(m,s)

  # Return the MLE
  return(theta)

}


#' Compute probability inverse transform values for Normal distribution
#'
#' @param x a numeric vector of length n
#'
#' @param theta a numeric vector of length two, containing MLE of parameters in Normal dist.
#'
#' @return a vector of length n containing the probability inverse transformed (PIT) values
#'
normalPIT = function(x, theta){

  # Extract MLE of mean and sd
  m <- theta[1]
  s <- theta[2]

  # Standardize values in the sample using MLE
  z <- (x - m) / s

  # Compute PIT values
  res <- pnorm(q = z, mean = 0, sd = 1)

  # Return the numeric vector
  return(res)
}


#' Compute Fisher information matrix by the negative expected value of Hessian matrix in Normal distribution.
#'
#' @param theta a numeric vector of length two, containing MLE of parameters in Normal dist
#'
#' @return Fisher information matrix for Normal distribution
#'
normalFisherByHessian = function(theta){

  # Extract the MLE of sd from theta vector
  sd.hat   <- theta[2]

  # Define a 2x2 matrix
  res      <- matrix(0, nrow = 2, ncol = 2)

  # Compute values
  res[1,1] <- 1/(sd.hat^2)
  res[1,2] <- res[2,1] <- 0
  res[2,2] <- 2/(sd.hat^2)

  # Return the Fisher information matrix
  return(res)

}
