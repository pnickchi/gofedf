#' Compute score function for Exponential distribution.
#'
#' @param x a numeric vector of length n
#'
#' @param theta a numeric vector of length two, containing MLE of parameters in Exponential dist.
#'
#' @return Score matrix with n rows and two columns.
#'
#' @noRd
#'
expScore = function(x, theta){

  # Compute score
  S <- (1/theta) - x

  # Convert the vector to a matrix
  S <- as.matrix(S)

  # Return the score
  return(S)
}



#' Calculate MLE of rate in Exponential dist.
#'
#' @param x a numeric vector of length n
#'
#' @return numeric value for the MLE of rate
#'
#' @noRd
#'
expMLE = function(x){
  theta <- 1 / mean(x)
  return(theta)
}



#' Compute probability inverse transform values for Exponential dist.
#'
#' @param x a numeric vector of length n
#'
#' @param theta a numeric vector of length one, containing MLE rate in Exponential dist.
#'
#' @return a vector of length n containing the probability inverse transformed (PIT) values
#'
#' @noRd
#'
expPIT = function(x, theta){

  res <- pexp(q = x, rate = theta)
  return(res)

}


