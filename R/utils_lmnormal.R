#' Compute score function for linear models.
#'
#' @param x a matrix with n rows and p columns containing the explanatory variables.
#'
#' @param y a numeric vector of length n containing the response variable.
#'
#' @param theta a numeric vector of length (p+1), containing MLE of parameters in a linear model.
#'
#' @return Score matrix with n rows and (p+1) columns.
#'
lmScore = function(x, y, theta){

  # Extract the MLE of coefficients in the linear model
  betahat <- theta$coef
  sigma2  <- theta$sigma2
  sigma   <- sqrt(sigma2)

  # Get the number of rows and columns in x matrix
  n       <- nrow(x)
  p       <- ncol(x)

  # Compute the fitted values
  yhat    <- x %*% betahat

  # Define a matrix for score
  S       <- matrix(NA, nrow = n, ncol = p + 1)

  # Compute the score matrix
  S[,1:p] <- ( x * as.numeric(y-yhat) ) / sigma2
  S[,p+1] <- (-1/sigma) + (y - yhat)^2 / sigma^3

  # Return the score matrix
  return(S)

}


#' Compute maximum likelihood estimates for linear models
#'
#' @param x a matrix with n rows and p columns containing the explanatory variables.
#'
#' @param y a numeric vector of length n containing the response variable.
#'
#' @return a numeric vector of estimates.
#'
lmMLE = function(x, y){

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


#' Compute probability inverse transform values for linear models.
#'
#' @param x a matrix with n rows and p columns containing the explanatory variables.
#'
#' @param y a numeric vector of length n containing the response variable.
#'
#' @param theta a numeric vector of length (p+1), containing MLE of parameters in a linear model.
#'
#' @return a vector of length n containing the probability inverse transformed (PIT) values
#'
lmPIT = function(x, y, theta){

  # Compute PIT values with pnorm function
  pit <- pnorm( (y - x %*% theta$coef ) / sqrt(theta$sigma2), mean = 0, sd = 1)
  pit <- as.numeric(pit)

  # Return pit values
  return(pit)

}


#' Compute Fisher information matrix in the case of linear model with Normal residuals.
#'
#' @param x a matrix with n rows and p columns containing the explanatory variables.
#'
#' @param y a numeric vector of length n containing the response variable.
#'
#' @param theta a numeric vector of length (p+1), containing MLE of parameters in a linear model.
#'
#' @return Fisher information matrix for linear models.
#'
lmFisherByHessian = function(x, y, theta){

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
