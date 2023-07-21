#' Apply Normal distribution to sample and compute required components for the test.
#'
#' @description Compute Maximum likelihood estimates of the parameters in Gamma distribution, Score function evaluated at the sample,
#' and probability inverse transformed (PIT) values of sample.
#'
#' @param x a numeric vector.
#'
#' @return a list with three elements:
#'
#' @noRd
applyNormal = function(x){

  # Compute MLE of parameters for Normal distribution
  n  <- length(x)
  m  <- mean(x)
  s  <- sqrt( (n-1) * var(x) / n )
  par <- c(m, s)

  # Compute score function
  S1 <- (x - m) / (s^2)
  S2 <- (x - m)^2/s^3 - rep(1/s,n)
  S  <- cbind(S1,S2)

  # Compute probability inverse transform of data
  pit <- pnorm( (x - m) / s, mean = 0, sd = 1)

  # Define the list to return
  res <- list(Score = S, pit = pit, par = par)
  return(res)
}
