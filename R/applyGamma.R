applyGamma = function(x, use.rate){

  # Calculate MLE of parameters in Gamma distribution
  par    <- estimateGamma(x, ur = use.rate)
  alpha  <- par[1]
  lambda <- par[2]

  # Compute score function for sample
  S1     <- log(lambda) - digamma(alpha) + log(x)
  S2     <- (alpha / lambda) - x
  S      <- cbind(S1,S2)

  # Calculate the probability inverse transfer of sample
  pit <- pgamma(q = x, shape = alpha, rate = lambda)

  # Define the list to return
  res <- list(Score = S, pit = pit, par = par)
  return(res)

}
