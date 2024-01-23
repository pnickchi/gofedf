#' Apply Normal distribution to sample and compute required components for the test.
#'
#' @description Compute maximum likelihood estimates (MLE) of the parameters in Normal distribution, score function evaluated at the MLE, and probability integral transformed (PIT) values of sample.
#'
#' @param x a numeric vector.
#'
#' @return a list with three elements:
#'
#' @noRd
applyNormal = function(x){

  # Call normalMLE() function to compute MLE of parameters for Normal distribution
  par <- normalMLE(x)

  # Call normalScore() function to compute score function
  S   <- normalScore(x = x, theta = par)

  # Call normalPIT() function to compute probability integral transformed values
  pit <- normalPIT(x = x, theta = par)

  # Define a list to return
  res <- list(Score = S, pit = pit, par = par)
  return(res)
}
