getFisherInformationForGLMGamma = function(par, x){

  n           <- nrow(x)
  rindx       <- ncol(x) + 1
  FI          <- matrix(NA, nrow = ncol(x)+2, ncol = ncol(x)+2)
  FI[1,1]     <- trigamma(par) - (1/par)
  FI[2:rindx,2:rindx] <- (-par/n) * ( t(X) %*% X )
  FI[1,2:rindx]   <- FI[2:rindx,1] <- 0

  return(FI)

}
