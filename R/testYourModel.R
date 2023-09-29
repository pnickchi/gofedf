#' Apply the Goodness of Fit Test Based on Empirical Distribution Function to Any Likelihood Model.
#'
#' @description This function applies the goodness-of-fit test based on empirical distribution function.
#' It requires certain inputs depending on whether the model involves parameter estimation or not.
#' If the model is known and there is no parameter estimation, the function requires the sample as a vector and
#' the probability transformed (or pit) values of the sample. This ought to be a vector as well. If there is
#' parameter estimation in the model, the function additionally requires the score as a matrix with n
#' rows and p columns, where n is the sample size and p is the number of estimated parameters.
#' The function checks if the score is zero at the estimated parameter (which is assumed to be the maximum
#' likelihood estimate).
#'
#' @param x a non-empty numeric vector of sample data.
#'
#' @param pit The probability transformed (or pit) values of the sample which ought to be a numeric vector with
#' the same size as x.
#'
#' @param score The default value is null and refers to no parameter estimation case. If there is parameter estimation,
#' the score matrix must be a matrix with n rows and p columns, where n is the sample size and p is the number of
#' estimated parameters.
#'
#' @param ngrid the number of equally spaced points to discretize the (0,1) interval for computing the covariance function.
#'
#' @param gridpit logical. If \code{TRUE} (the default value), the parameter ngrid is ignored and (0,1) interval is divided
#' based on probability inverse transformed values obtained from the sample. If \code{FALSE}, the interval is divided into ngrid
#' equally spaced points for computing the covariance function.
#'
#' @param precision The theory behind goodness-of-fit test based on empirical distribution function (edf) works well
#' if the MLE is indeed the root of derivative of log likelihood function. A precision of 1e-9 (default value) is used
#' to check this. A warning message is generated if the score evaluated at MLE is not close enough to zero.
#'
#' @param method a character string indicating which goodness-of-fit statistic is to be computed. The default value is
#' 'cvm' for the Cramer-von-Mises statistic. Other options include 'ad' for the Anderson-Darling statistic, and 'both'
#' to compute both cvm and ad.
#'
#'
#' @return A list of two containing the following components:
#' - Statistic: the value of goodness-of-fit statistic.
#' - p-value: the approximate p-value for the goodness-of-fit test based on empirical distribution function.
#' if method = 'cvm' or method = 'ad', it returns a numeric value for the statistic and p-value. If method = 'both', it
#' returns a numeric vector with two elements and one for each statistic.
#'
#' @export
#'
#' @import stats
#' @import statmod
#'
#' @examples
#' # Example: Inverse Gaussian (IG) distribution with weights
#'
#' # Set the seed to reproduce example.
#' set.seed(123)
#'
#' # Set the sample size
#' n <- 50
#'
#' # Assign weights
#' weights <- rep(1.5,n)
#'
#' # Set mean and shape parameters for IG distribution.
#' mio        <- 2
#' lambda     <- 2
#'
#' # Generate a random sample from IG distribution with weighted shape.
#' sim_data <- statmod::rinvgauss(n, mean = mio, shape = lambda * weights)
#'
#' # Compute MLE of parameters, score matrix, and pit values.
#' theta_hat    <- inversegaussianMLE(obs = sim_data,   w = weights)
#' ScoreMatrix  <- inversegaussianScore(obs = sim_data, w = weights, mle = theta_hat)
#' pitvalues    <- inversegaussianPIT(obs = sim_data ,  w = weights, mle = theta_hat)
#'
#' # Apply the goodness-of-fit test.
#' testYourModel(x = sim_data, pit = pitvalues, score = ScoreMatrix)
#'
testYourModel = function(x, pit, score = NULL, ngrid = length(x), gridpit = TRUE, precision = 1e-9, method = 'cvm'){

  if( !is.double(x) ){
    stop('x values must be numeric.')
  }

  if( !is.vector(x) ){
    stop('x must be a vector of numeric values.')
  }

  if( length(x) < 2 ){
    stop('The number of observations in x must be greater than or equal two.')
  }

  if( any(pit < 0) | any(pit > 1) ){
    stop('pit values must be between zero and one.')
  }

  #
  # Case 1: if score is null then it means there was no parameter estimation in the model.
  # The model is fully specified and the covariance function of the  stochastic process,
  # $\hat W_{n}(u)$, is simply min(s,t)-st for 0 <= s,t <=1. The Eigen values are:
  # (i) for Cramer-von-Mises:  1 / ( pi^2 * j^2 ); j=1,2,3...
  # (ii) for Anderson-Darling: 1 / ( j * (j+1) ); j=1,2,3,...
  #
  if( is.null(score) ){

    if ( method == 'cvm') {

      # Calculate Cramer-von-Mises statistic
      cvm        <- getCvMStatistic(pit)
      names(cvm) <- 'Cramer-von-Mises Statistic'

      # Compute the 100 leading Eigenvalues of the covariance function
      j  <- 1:100
      ev <- 1 / ( pi^2 * j^2 )

      # Calculate pvalue
      pvalue  <- getpvalue(u = cvm, eigen = ev)

      # Prepare a list to return statistic and pvalue
      res     <- list(Statistic = cvm, pvalue = pvalue)

      return(res)

    } else if ( method == 'ad') {

      # Calculate Anderson-Darling statistic
      AD        <- getADStatistic(pit)
      names(AD) <- 'Anderson-Darling Statistic'

      # Compute the 100 leading Eigenvalues of the covariance function
      j <- 1:100
      ev <- 1 / ( j * (j+1) )

      # Calculate pvalue
      pvalue  <- getpvalue(u = AD, eigen = ev)

      # Prepare a list to return statistic and pvalue
      res     <- list(Statistic = AD, pvalue = pvalue)
      return(res)

    } else {

      # Calculate both cvm and ad statisitcs

      # 1. Do cvm calculation
      cvm        <- getCvMStatistic(pit)
      names(cvm) <- 'Cramer-von-Mises Statistic'

      # Compute the 100 leading Eigenvalues of the covariance function
      j <- 1:100
      ev <- 1 / ( pi^2 * j^2 )

      # Calculate pvalue
      cvm.pvalue  <- getpvalue(u = cvm, eigen = ev)
      names(cvm.pvalue) <- 'pvalue for Cramer-von-Mises test'

      # 2. Do AD calculations
      ad      <- getADStatistic(pit)
      names(ad) <- 'Anderson-Darling Statistic'

      # Compute the 100 leading Eigenvalues of the covariance function
      j <- 1:100
      ev <- 1 / ( j * (j+1) )

      # Calculate pvalue
      ad.pvalue  <- getpvalue(u = ad, eigen = ev)
      names(ad.pvalue) <- 'Anderson-Darling test'

      # Prepare a list to return both statistics and their approximate pvalue
      res     <- list(Statistics = c(cvm, ad), pvalue = c(cvm.pvalue, ad.pvalue) )
      return(res)
    }

  }

  #
  # Case 2: if score is not null then it means there is parameter estimation in the model.
  # The score should be a matrix with n rows (number of observations in x) and p columns
  # (the number of estimated parameters in the model).
  # The covariance function of the stochastic process, $\hat W_{n}(u)$, needs to be calculated
  # based on the GoF based on EDF.
  #
  if( !is.null(score) ){

    # Check if the score is a matrix
    if( !is.matrix(score) ){
      stop('score must be a matrix.')
    }

    # Get the length of x.
    n    <- length(x)

    # Check if score dimension is correct
    if( nrow(score) != n ){
      stop('The number of rows in score matrix do not match the sample size in x.')
    }

    # Check if score is zero at MLE
    if( any(colSums(score) > precision) ){
      warning( paste0('Score matrix is not zero at MLE. precision of ', precision, ' was used') )
    }

    # Calculate Fisher information matrix by computing the variance of score from the sample.
    fisher  <- (n-1)*var(score)/n


    # Calculate Cramer-von-Mises, Anderson-Darling statistics or both
    if ( method == 'cvm' ) {

      # Calculate Cramer-von-Mises statistic
      cvm        <- getCvMStatistic(pit)
      names(cvm) <- 'Cramer-von-Mises Statistic'

      # Get Eigen values
      if( gridpit ){
        ev    <- getEigenValues(S = score, FI = fisher, pit = pit, me = 'cvm')
      }else{
        ev    <- getEigenValues_manualGrid(S = score, FI = fisher, pit = pit, M = ngrid, me = 'cvm')
      }

      # Calculate pvalue
      pvalue  <- getpvalue(u = cvm, eigen = ev)

      # Prepare a list to return statistic and pvalue
      res     <- list(Statistic = cvm, pvalue = pvalue)

      return(res)

    } else if ( method == 'ad') {

      # Calculate Anderson-Darling statistic
      AD      <- getADStatistic(pit)
      names(AD) <- 'Anderson-Darling Statistic'

      # Get Eigen values
      if( gridpit ){
        ev    <- getEigenValues(S = score, FI = fisher, pit = pit, me = 'ad')
      }else{
        ev    <- getEigenValues_manualGrid(S = score, FI = fisher, pit = pit, M = ngrid, me = 'ad')
      }

      # Calculate pvalue
      pvalue  <- getpvalue(u = AD, eigen = ev)

      # Prepare a list to return statistic and pvalue
      res     <- list(Statistic = ad, pvalue = pvalue)

      return(res)

    } else {

      # Calculate both cvm and ad statisitcs

      # 1. Do cvm calculation
      cvm        <- getCvMStatistic(pit)
      names(cvm) <- 'Cramer-von-Mises Statistic'

      # Get Eigen values
      if( gridpit ){
        ev    <- getEigenValues(S = score, FI = fisher, pit = pit, me = 'cvm')
      }else{
        ev    <- getEigenValues_manualGrid(S = score, FI = fisher, pit = pit, M = ngrid, me = 'cvm')
      }

      # Calculate pvalue
      cvm.pvalue  <- getpvalue(u = cvm, eigen = ev)
      names(cvm.pvalue) <- 'pvalue for Cramer-von-Mises test'


      # 2. Do ad calculations
      AD        <- getADStatistic(pit)
      names(AD) <- 'Anderson-Darling Statistic'

      # Get Eigen values
      if( gridpit ){
        ev    <- getEigenValues(S = score, FI = fisher, pit = pit, me = 'ad')
      }else{
        ev    <- getEigenValues_manualGrid(S = score, FI = fisher, pit = pit, M = ngrid, me = 'ad')
      }

      # Calculate pvalue
      ad.pvalue  <- getpvalue(u = AD, eigen = ev)
      names(ad.pvalue) <- 'Anderson-Darling test'


      # Prepare a list to return both statistics and their approximate pvalue
      res     <- list(Statistics = c(cvm, AD), pvalue = c(cvm.pvalue, ad.pvalue) )
      return(res)

    }

  }

}
