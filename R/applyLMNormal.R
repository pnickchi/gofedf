applyLMNormal = function(x, y){

  # Calculate MLE of parameters
  par <- getMLEinLMNormal(x = x, y = y)

  # Get Score function
  Score <- getScoreinLMNormal(x, y, theta = par)

  # Get probability inverse transformation (pit)
  pit <- pnorm( (y - x %*% par$coef ) / sqrt(par$sigma2), mean = 0, sd = 1)
  pit <- as.numeric(pit)

  # Return the list
  return( list(par = par, Score = Score, pit = pit) )
}
