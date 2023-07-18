applyExponential = function(x){

  # Calculate MLE of parameter for Exponential distribution
  n      <- length(x)
  lambda <- 1 / mean(x)
  par    <- lambda

  # Compute score function for sample
  S1 <- (1/lambda) - x
  S  <- cbind(S1)

  # Calculate the probability inverse transform of sample
  pit <- pexp(q = x, rate = lambda, lower.tail = TRUE)

  # Define the list to return
  res <- list(Score = S, pit = pit, par = par)
  return(res)
}
