applyNormal                 = function(x){

  # Calculate MLE of parameters for normal distribution
  n  <- length(x)
  m  <- mean(x)
  s  <- sqrt( (n-1) * var(x) / n )
  par <- c(m, s)

  # Define score function matrix
  S1 <- (x - m) / (s^2)
  S2 <- (x - m)^2/s^3 - rep(1/s,n)
  S  <- cbind(S1,S2)

  # Calculate probability inverse transform of data
  pit <- pnorm( (x - m) / s, mean = 0, sd = 1)

  # Define the list to return
  res <- list(Score = S, pit = pit, par = par)
  return(res)
}
