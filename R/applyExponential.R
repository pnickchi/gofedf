#' Apply Exponential distribution to sample and compute required components for the test.
#'
#' @description Compute Maximum likelihood estimates of the parameters in Exponential distribution, Score function evaluated at the sample,
#' and probability inverse transformed (PIT) values of sample.
#'
#' @param x a numeric vector.
#'
#' @return a list with three elements.
#'
#' @noRd
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
