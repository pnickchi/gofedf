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
#' @noRd
applyGamma = function(x, use.rate){

  # Compute MLE of parameters in Gamma distribution
  par    <- estimateGamma(x, ur = use.rate)
  alpha  <- par[1]
  lambda <- par[2]

  # Compute score function for sample
  S1     <- log(lambda) - digamma(alpha) + log(x)
  S2     <- (alpha / lambda) - x
  S      <- cbind(S1,S2)

  # Compute the probability inverse transfer of sample
  pit <- pgamma(q = x, shape = alpha, rate = lambda)

  # Define the list to return
  res <- list(Score = S, pit = pit, par = par)
  return(res)

}
