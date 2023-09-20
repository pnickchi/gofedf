#' Apply linear model and extract required components for the test
#'
#' @param x is either a numeric vector or a design matrix. In the design matrix, rows indicate observations and columns
#' presents explanatory variables.
#'
#' @param y is a vector of numeric values with the same number of observations or number of rows as x.
#'
#' @return a list with three elements.
#'
applyLMNormal = function(x, y){

  # Compute MLE of parameters
  par <- lmMLE(x = x, y = y)

  # Compute score function for sample
  Score <- lmScore(x, y, theta = par)

  # Compute the probability inverse transfer (pit) values
  pit <- lmPIT(x = x, y = y, theta = par)

  # pit <- pnorm( (y - x %*% par$coef ) / sqrt(par$sigma2), mean = 0, sd = 1)
  # pit <- as.numeric(pit)

  # Return the list
  return( list(par = par, Score = Score, pit = pit) )
}
