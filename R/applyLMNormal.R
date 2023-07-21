#' Apply linear model to sample with the assumption of normal response and compute required components for the test.
#'
#' @param x a matrix of explanatory variables.
#' @param y a numeric vector of response values.
#'
#' @return a list with three elements.
#'
#' @noRd
applyLMNormal = function(x, y){

  # Compute MLE of parameters
  par <- getMLEinLMNormal(x = x, y = y)

  # Compute score function for sample
  Score <- getScoreinLMNormal(x, y, theta = par)

  # Compute the probability inverse transfer of sample
  pit <- pnorm( (y - x %*% par$coef ) / sqrt(par$sigma2), mean = 0, sd = 1)
  pit <- as.numeric(pit)

  # Return the list
  return( list(par = par, Score = Score, pit = pit) )
}
