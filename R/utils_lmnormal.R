#' Compute MLE estimates in the case of linear model.
#'
#' @param x a matrix with explanatory variables.
#'
#' @param y a numeric vector of response variables.
#'
#' @return a numeric vector of estimates.
#'
#' @noRd
getMLEinLMNormal = function(x, y){

  # Get sample size
  n          <- nrow(x)

  # Estimate MLE of coefficient
  coefhat    <- solve( t(x) %*% x ) %*% t(x) %*% y
  coefhat    <- as.vector(coefhat)

  # Estimate MLE of error terms
  sigma2hat  <- t(y - x %*% coefhat) %*% (y - x %*% coefhat)
  sigma2hat  <- as.numeric( sigma2hat / n )

  return( list(coef = coefhat, sigma2 = sigma2hat) )
}


#' Compute score function in the case of linear model.
#'
#' @param x a matrix with explanatory variables.
#'
#' @param y a numeric vector of response variables.
#'
#' @param theta a numeric vector with mle values.
#'
#' @return a score matrix, rows present samples and columns presents the variables.
#'
#' @noRd
getScoreinLMNormal = function(x, y, theta){

  # Extract the mle of coefficient
  betahat <- theta$coef
  sigma2  <- theta$sigma2
  sigma   <- sqrt(sigma2)

  # Get the number of rows and columns in x matrix
  n       <- nrow(x)
  p       <- ncol(x)

  # Find the fitted values
  yhat    <- x %*% betahat

  # Compute the score matrix
  S       <- matrix(NA, nrow = n, ncol = p + 1)
  S[,1:p] <- ( x* as.numeric(y-yhat) ) / sigma2
  S[,p+1] <- (-1/sigma) + (y - yhat)^2 / sigma^3

  # Return the score matrix
  return(S)

}


#' Compute Fisher information matrix in the case of linear model with Normal residuals.
#'
#' @param x a matrix with explanatory variables.
#'
#' @param y a numeric vector of response variables.
#'
#' @param theta a numeric vector with mle values.
#'
#' @return Fisher information matrix.
#'
#' @noRd
#'
getObservedHessMatrixinLMNormal = function(x, y, theta){

  sigma2   <- theta$sigma2

  p <- ncol(x)
  FI <- matrix(NA, nrow = p+1, ncol = p+1)

  FI[1:p,1:p] <- t(x)%*%x / sigma2

  FI[1:p, p+1] <- 0
  FI[p+1,1:p]  <- 0

  FI[1:p,p+1] <- 0
  FI[p+1,p+1] <- (2*n/sigma2)

  return(FI)
}
