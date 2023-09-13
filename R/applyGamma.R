#' Apply Gamma distribution to sample and compute required components for the test.
#'
#' @description Compute Maximum likelihood estimates of the parameters in Gamma distribution, Score function evaluated at the sample,
#' and probability inverse transformed (PIT) values of sample.
#'
#' @param x a numeric vector.
#'
#' @param use.rate logical. If \code{TRUE} the rate parameter is returned while estimating MLE. Otherwise the scale is returned.
#'
#' @return a list with three elements.
#'
applyGamma = function(x, use.rate){

  # Compute MLE of parameters in Gamma distribution
  par    <- gammaMLE(x, ur = use.rate)

  # Compute score function for sample
  S <- gammaScore(x = x, theta = par)

  # Compute the probability inverse transfer of sample
  pit <- gammaPIT(x = x, theta = par)

  # Define the list to return
  res <- list(Score = S, pit = pit, par = par)
  return(res)

}
