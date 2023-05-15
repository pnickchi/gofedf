#' Compute the maximum likelihood estimate of parameters in Inverse Gaussian distribution with weighted observations.
#'
#' @param obs a numeric vector of sample observations.
#' @param ... a list of additional parameters to define the likelihood. In this function, weight is being passed.
#'
#' @export
#'
#' @return The function compute the MLE of parameters in Inverse Gaussian distribution and returns a vector of
#' estimates. The first and second elements of the vector are MLE of the mean and MLE of shape, respectively.
#'
IG_mlefunc = function(obs, ...){

  args <- list(...)
  w <- args$w

  if( any(w<0) ){
    stop('The weights must be positive.')
  }

  if( any(obs<0) ){
    stop('Negative values are not allowed in Inverse Guassian distribution.')
  }

  n      <- length(obs)
  mio    <- sum(w * obs) / sum(w)
  lambda <- n / ( sum( (w/obs) - (w/mio) ) )
  return(c(mio,lambda))

}


#' Compute the score function of the Inverse Gaussian distribution based on a sample.
#'
#' @param obs a numeric vector of sample observations.
#' @param ... a list of additional parameters to define the likelihood. In this function, weight and mle are being passed.
#'
#' @export
#'
#' @return The score matrix with n rows (number of sample observations) and 2 columns (mean and shape).
IG_scorefunc = function(obs, ...){

  args <- list(...)
  mle <- args$mle
  w   <- args$w

  # Get MLE estimates of mio and lambda
  mio    <- mle[1]
  lambda <- mle[2]

  # Assign x as weights and y as sample data.
  x <- w
  y <- obs


  # First column of score matrix (mean)
  S1 <- ( (lambda * x * y)/mio^3 ) - ( (lambda * x)/mio^2 )

  # Second column of score matrix (shape)
  S2 <- (1/(2*lambda)) + (x/mio) - ( (x*y) / (2*mio^2) ) - (x / (2*y))

  # Create score matrix
  S <- cbind(S1,S2)

  # return score
  return(score = S)
}



#' Compute the probability transformed values for a sample from Inverse Gaussian distribution.
#'
#' @param obs A numeric vector of sample observations.
#' @param ... A list of additional parameters to define the likelihood. In this function, MLE is being passed.
#'
#' @export
#'
#' @import statmod
#'
#' @return A numeric vector of probability transformed values of sample observations.
IG_pitfunc = function(obs, ...){

  # Extract the MLE
  args <- list(...)
  mle <- args$mle

  # Set observations as y
  y <- obs

  # Compute pit values with the estimated MLE
  pit <- statmod::pinvgauss(q = y, mean = mle[1], shape = mle[2])

  # return pit values
  return(pit)

}