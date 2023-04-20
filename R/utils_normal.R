#' Score function for Normal dist
#'
#' @param x vector of data
#' @param theta mle
#'
#' @return Score matrix
#' @export
#'
normalScore = function(x, theta){

  n <- length(x)

  m <- theta[1]
  s <- theta[2]

  S1 <- (x - m) / (s^2)
  S2 <- (x - m)^2/s^3 - rep(1/s, n)
  S  <- cbind(S1,S2)

  return(S)
}


#' Get MLE estimate for Normal
#'
#' @param x vector of data
#'
#' @return MLE estimates of mean and sd
#' @export
#'
normalMLE = function(x){

  n <- length(x)
  m <- mean(x)
  s <- sqrt( (n-1)* var(x) / n)
  theta <- c(m,s)
  return(theta)

}


#' Calculate probability inverse transform for Normal distribution
#'
#' @param x vector of data
#' @param theta mle
#'
#' @return vector of probability inverse transformed values
#' @export
#'
normalPIT = function(x, theta){
  z <- (x - theta[1]) / theta[2]
  res <- pnorm(q = z, mean = 0, sd = 1)
  return(res)
}



#' Calculate observed Hessian matrix for Normal distribution
#'
#' @param par a vector of mle parameters
#'
#' @return Observed Hessian matrix for Normal distribution
#' @export
#'
observedHessianMatrixNormal = function(par){

  sd.hat   <- par[2]
  res      <- matrix(0, nrow = 2, ncol = 2)
  res[1,1] <- 1/(sd.hat^2)
  res[1,2] <- res[2,1] <- 0
  res[2,2] <- 2/(sd.hat^2)
  return(res)

}
