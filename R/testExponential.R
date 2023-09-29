#' Apply Goodness of Fit Test for Exponential Distribution
#'
#' @description Performs the goodness-of-fit test based on empirical distribution function to check if an i.i.d sample
#' follows an Exponential distribution.
#'
#' @param x a non-empty numeric vector of sample data.
#'
#' @param ngrid the number of equally spaced points to discretize the (0,1) interval for computing the covariance function.
#'
#' @param gridpit logical. If \code{TRUE} (the default value), the parameter ngrid is ignored and (0,1) interval is divided
#' based on probability inverse transformed values obtained from the sample. If \code{FALSE}, the interval is divided into ngrid
#' equally spaced points for computing the covariance function.
#'
#' @param hessian logical. If \code{TRUE} the Fisher information matrix is estimated by the observed Hessian Matrix based on
#' the sample. If \code{FALSE} (the default value) the Fisher information matrix is estimated by the variance of the
#' observed score matrix.
#'
#' @param method a character string indicating which goodness-of-fit statistic is to be computed. The default value is
#' 'cvm' for the Cramer-von-Mises statistic. Other options include 'ad' for the Anderson-Darling statistic, and 'both'
#' to compute both cvm and ad.
#'
#' @return A list of two containing the following components:
#' - Statistic: the value of goodness-of-fit statistic.
#' - p-value: the approximate p-value for the goodness-of-fit test based on empirical distribution function.
#' if method = 'cvm' or method = 'ad', it returns a numeric value for the statistic and p-value. If method = 'both', it
#' returns a numeric vector with two elements and one for each statistic.
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' n <- 50
#' sim_data <- rexp(n, rate = 2)
#' testExponential(x = sim_data)
testExponential = function(x, ngrid = length(x), gridpit = FALSE, hessian = FALSE, method = 'cvm'){

  # Check if any observation is negative
  if( any(x < 0) ){
    stop('x values must be positive for Exponential distribution.')
  }

  # Get the number of sample
  n       <- length(x)

  # Apply exponential distribution
  temp    <- applyExponential(x)

  # Extract score function, pit values, and MLE estimate
  Score   <- temp$Score
  pit     <- temp$pit
  par     <- temp$par

  # Compute Fisher information matrix
  if( hessian ){
    fisher <- matrix( n / (par[1])^2, nrow = 1, ncol = 1 )
  }else{
    fisher <- (n-1)*var(Score)/n
  }

  if( method == 'cvm'){

    # Compute Cramer-von-Mises statistic
    cvm      <- getCvMStatistic(pit)

    # Compute Eigen values
    if( gridpit ){
      ev    <- getEigenValues(S = Score, FI = fisher, pit, me = 'cvm')
    }else{
      ev    <- getEigenValues_manualGrid(S = Score, FI = fisher, pit, M = ngrid, me = 'cvm')
    }

    # Compute p-value
    pvalue  <- getpvalue(u = cvm, eigen = ev)

    res     <- list(Statistic = cvm, pvalue = pvalue)
    return(res)

  } else if ( method == 'ad') {

    # Compute Anderson-Darling statistic
    AD      <- getADStatistic(pit)

    # Compute Eigen values
    if( gridpit ){
      ev    <- getEigenValues(S = Score, FI = fisher, pit, me = 'ad')
    }else{
      ev    <- getEigenValues_manualGrid(S = Score, FI = fisher, pit, M = ngrid, me = 'ad')
    }

    # Compute pvalue
    pvalue  <- getpvalue(u = AD, eigen = ev)
    res     <- list(Statistic = AD, pvalue = pvalue)

    return(res)

  } else {

    # Calculate both cvm and ad statistics

    # Compute Cramer-von-Mises statistic
    cvm        <- getCvMStatistic(pit)
    names(cvm) <- 'Cramer-von-Mises Statistic'

    # Compute Eigen values
    if( gridpit ){
      ev    <- getEigenValues(S = Score, FI = fisher, pit = pit, me = 'cvm')
    }else{
      ev    <- getEigenValues_manualGrid(S = Score, FI = fisher, pit = pit, M = ngrid, me = 'cvm')
    }

    # Compute pvalue
    cvm.pvalue  <- getpvalue(u = cvm, eigen = ev)
    names(cvm.pvalue) <- 'pvalue for Cramer-von-Mises test'


    # Compute Anderson-Darling statistic
    ad      <- getADStatistic(pit)
    names(ad) <- 'Anderson-Darling Statistic'

    # Compute Eigen values
    if( gridpit ){
      ev    <- getEigenValues(S = Score, FI = fisher, pit = pit, me = 'ad')
    }else{
      ev    <- getEigenValues_manualGrid(S = Score, FI = fisher, pit = pit, M = ngrid, me = 'ad')
    }

    # Compute pvalue
    ad.pvalue  <- getpvalue(u = ad, eigen = ev)
    names(ad.pvalue) <- 'Anderson-Darling test'

    # Prepare a list to return both statistics and their approximate pvalue
    res     <- list(Statistics = c(cvm, ad), pvalue = c(cvm.pvalue, ad.pvalue) )
    return(res)

  }

}
