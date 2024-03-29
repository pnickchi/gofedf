#' Apply Goodness of Fit Test for Gamma Distribution
#'
#' @description Performs the goodness-of-fit test based on empirical distribution function to check if an i.i.d sample
#' follows a Gamma distribution.
#'
#' @param x a non-empty numeric vector of sample data.
#'
#' @param ngrid the number of equally spaced points to discretize the (0,1) interval for computing the covariance function.
#'
#' @param gridpit logical. If \code{TRUE} (the default value), the parameter ngrid is ignored and (0,1) interval is divided
#' based on probability integral transformed values obtained from the sample. If \code{FALSE}, the interval is divided into ngrid
#' equally spaced points for computing the covariance function.
#'
#' @param hessian logical. If \code{TRUE} the Fisher information matrix is estimated by the observed Hessian Matrix based on
#' the sample. If \code{FALSE} (the default value) the Fisher information matrix is estimated by the variance of the
#' observed score matrix.
#'
#' @param rate logical. If \code{TRUE} (the default value), the rate is estimated in Gamma distribution. If \code{FALSE}, scale
#' is estimated. See \code{\link{GammaDist}} for more details.
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
#' sim_data <- rgamma(n = 50, shape = 3)
#' testGamma(x = sim_data)
#' sim_data <- runif(n = 50)
#' testGamma(x = sim_data)
testGamma = function(x, ngrid = length(x), gridpit = FALSE, hessian = FALSE, rate = TRUE, method = 'cvm'){


  if( !is.numeric(x) | !is.vector(x) ){
    stop('x must be a numeric vector.')
  }

  # Make sure all observations are positive
  if( any(x < 0) ){
    stop('x values must be positive for Gamma distribution')
  }

  if( !(ngrid > 0) ){
    stop('ngrid must be a positive number.')
  }

  if( !(ngrid %% 1 == 0) ){
    stop('ngrid must be an integer number.')
  }

  if( !is.logical(gridpit) ){
    stop('gridpit must be either TRUE or FALSE.')
  }

  if( !is.logical(hessian) ){
    stop('hessian must be either TRUE or FALSE.')
  }

  if( !is.vector(method) | length(method) > 1){
    stop('method must be a character string with length one.')
  }

  if( !(method %in% c('cvm','ad','both')) ){
    stop('ethod must be either cvm, ad, or both.')
  }

  # Get sample size
  n      <- length(x)

  # Apply Gamma distribution
  temp   <- applyGamma(x, use.rate = rate)

  # Extract score function, pit values, and MLE
  Score  <- temp$Score
  pit    <- temp$pit
  par    <- temp$par

  # Compute Fisher information matrix
  if( hessian ){
    fisher <- gammaFisherByHessian(par)
  }else{
    fisher <- (n-1)*var(Score)/n
  }


  if( method == 'cvm'){

    # Compute Eigen values
    if( gridpit ){
      ev    <- getEigenValues(S = Score, FI = fisher, pit, me = 'cvm')
    }else{
      ev    <- getEigenValues_manualGrid(S = Score, FI = fisher, pit, M = ngrid, me = 'cvm')
    }

    # Compute Cramer-von-Mises statistic
    cvm     <- getCvMStatistic(pit)

    # Compute pvalue
    pvalue  <- getpvalue(u = cvm, eigen = ev)

    res     <- list(Statistic = cvm, pvalue = pvalue)
    return(res)

  } else if ( method == 'ad') {


    # Compute Eigen values
    if( gridpit ){
      ev    <- getEigenValues(S = Score, FI = fisher, pit, me = 'ad')
    }else{
      ev    <- getEigenValues_manualGrid(S = Score, FI = fisher, pit, M = ngrid, me = 'ad')
    }

    # Compute Anderson-Darling statistic
    AD      <- getADStatistic(pit)

    # Compute pvalue
    pvalue  <- getpvalue(u = AD, eigen = ev)

    res     <- list(Statistic = AD, pvalue = pvalue)
    return(res)

  } else {

    # Calculate both cvm and ad statisitcs

    # 1. Do cvm calculation
    cvm        <- getCvMStatistic(pit)
    names(cvm) <- 'Cramer-von-Mises Statistic'

    # Get Eigen values
    if( gridpit ){
      ev    <- getEigenValues(S = Score, FI = fisher, pit = pit, me = 'cvm')
    }else{
      ev    <- getEigenValues_manualGrid(S = Score, FI = fisher, pit = pit, M = ngrid, me = 'cvm')
    }

    # Calculate pvalue
    cvm.pvalue  <- getpvalue(u = cvm, eigen = ev)
    names(cvm.pvalue) <- 'pvalue for Cramer-von-Mises test'


    # 2. Do ad calculations
    ad      <- getADStatistic(pit)
    names(ad) <- 'Anderson-Darling Statistic'

    # Get Eigen values
    if( gridpit ){
      ev    <- getEigenValues(S = Score, FI = fisher, pit = pit, me = 'ad')
    }else{
      ev    <- getEigenValues_manualGrid(S = Score, FI = fisher, pit = pit, M = ngrid, me = 'ad')
    }

    # Calculate pvalue
    ad.pvalue  <- getpvalue(u = ad, eigen = ev)
    names(ad.pvalue) <- 'Anderson-Darling test'

    # Prepare a list to return both statistics and their approximate pvalue
    res     <- list(Statistics = c(cvm, ad), pvalue = c(cvm.pvalue, ad.pvalue) )
    return(res)

  }

}
