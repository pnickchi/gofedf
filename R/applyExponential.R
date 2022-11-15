applyExponential = function(x){

  # Calculate MLE of parameter for Exp distribution
  n      <- length(x)
  lambda <- 1 / mean(x)
  par    <- lambda

  # Define score function matrix
  S1 <- (1/lambda) - x
  S  <- cbind(S1)

  # Calculate probability inverse transform of data
  pit <- pexp(q = x, rate = lambda, lower.tail = TRUE)

  # Define the list to return
  res <- list(Score = S, pit = pit, par = par)
  return(res)
}
