observedHessianMatrixGamma = function(par){

  alpha.hat  <- par[1]
  lambda.hat <- par[2]

  res      <- matrix(0, nrow = 2, ncol = 2)
  res[1,1] <- trigamma(alpha.hat)
  res[1,2] <- res[2,1] <- -1/lambda.hat
  res[2,2] <- alpha.hat / lambda.hat^2
  return(res)

}

estimateGamma <- function(x, ur){

  n <- length(x)
  m <- mean(x)  #m1
  s <- var(x)   #m2
  b <- s / m

  a <- m / b
  mlog <- mean(log(x))

  logm <- log(m)
  aold <- a
  anew <- aold -(log(aold)-logm + mlog -digamma(aold))/(1/aold-trigamma(aold))
  bnew=m/anew

  if( anew < 0) anew <- aold/2

  while ( abs(anew-aold) > 1e-7){
    aold <- anew
    old.score = (log(aold)-log(m)+ mlog -digamma(aold))
    old.score.derivative = 1/aold-trigamma(aold)
    anew <- aold - old.score/old.score.derivative
    if( anew < 0) anew <- aold/2
  }

  beta  <- m/anew
  alpha <- anew

  if(ur){
    return(c(alpha,1/beta))
  }
  return(c(alpha, beta))
}



# F(x) function for Gamma dist
GammaFx = function(x, theta){

  if( x > 0){
    res <- pgamma(q = x, shape = theta[1], rate = theta[2])
  }else{
    res <- 0
  }

  return(res)
}

# Score function for Gamma dist
GammaScore = function(x, theta){

  alpha  <- theta[1]
  lambda <- theta[2]

  S1 <- log(lambda) - digamma(alpha) + log(x)
  S2 <- (alpha / lambda) - x
  S  <- cbind(S1,S2)

  return(S)
}
