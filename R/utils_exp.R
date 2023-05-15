#' The Probability Transform Inverse Function for Exponential Distribution
#'
#' @param x a vector data
#' @param theta MLE estimate of data
#'
#' @return A vector of transformed values
#'
expFx = function(x, theta){

  if( all(x > 0) ){
    res <- 1 - exp(-theta*x)
  }else{
    res <- 0
  }

  return(res)

}

#' Calculate Score Function for Exponential Distribution
#'
#' @param x a vector data
#' @param theta MLE estimate of data
#'
#' @return A matrix of nx1
#'
expScore = function(x, theta){
  S <- (1/theta) - x
  S <- as.matrix(S)
  return(S)
}

#' Calculate MLE of Parameter in Exponential Distribution
#'
#' @param x a vector data
#'
#' @return Numeric
#'
expMLE = function(x){
  theta <- 1 / mean(x)
  return(theta)
}
