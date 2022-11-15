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

getScoreinLMNormal = function(x, y, theta){

  betahat <- theta$coef
  sigma2  <- theta$sigma2
  s       <- sqrt(sigma2)

  n <- nrow(x)
  p <- ncol(x)
  yhat <- x %*% betahat
  S <- matrix(NA, nrow = n, ncol = p + 1)
  S[,1:p] <- ( x* as.numeric(y-yhat) ) / sigma2
  S[,p+1] <- (-1/s) + (y - yhat)^2 / s^3

  return(S)

}

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
