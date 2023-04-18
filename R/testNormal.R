#' Apply Goodness of Fit Test for Normal Distribution
#'
#' @param x A numeric vector of data values with length n.
#' @param ngrid The number of equally spaced points to discritize the (0,1) interval to estimate the covariance of the process.
#' @param gridpit A Boolean indicator. If TRUE, ngrid is ignored and (0,1) interval is divided based on probability inverse transformed
#'  values. If FALSE (default value), (0,1) is divided into ngrid equally spaced points to estimate process.
#' @param hessian A Boolean indicator to control if Hessian matrix should be used in estimation of Fisher information matrix.
#' If False (the default value), the Fisher information matrix is estimated by the observed score function.
#' If TRUE, the Fisher information matrix is estimated by the observed Hessian Matrix.
#'
#' @param method a character string to indicate which statistics to calculate.
#' The possible values are 'cvm' for Cramer-von Mises and 'ad' for Anderson-Darling.
#' The default value is 'cvm'.
#'
#' @return A list of two.
#' - Statistic: The Cramer-von-Mises statistic.
#' - pvalue: The approximated pvalue for the GoF test based on EDF.
#' @export
#'
#' @examples
#' set.seed(123)
#' sim_data <- rnorm(n = 100)
#' testNormal(x = sim_data)
#' sim_data <- rgamma(100, shape = 3)
#' testNormal(x = sim_data)
testNormal = function(x, ngrid = length(x), gridpit = FALSE, hessian = FALSE, method = 'cvm'){

  if( !is.numeric(x) | !is.vector(x) ){
     stop('x must be a numeric vector.')
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
     stop('method must be either cvm, ad, or both.')
  }

  n       <- length(x)
  temp    <- applyNormal(x)
  Score   <- temp$Score
  pit     <- temp$pit
  par     <- temp$par

  if( hessian ){
    fisher <- observedHessianMatrixNormal(par)
  }else{
    fisher <- (n-1)*var(Score)/n
  }


  if( method == 'cvm'){

    cvm      <- getCvMStatistic(pit)

    if( gridpit ){
      ev    <- getEigenValues(S = Score, FI = fisher, pit, me = 'cvm')
    }else{
      ev    <- getEigenValues_manualGrid(S = Score, FI = fisher, pit, M = ngrid, me = 'cvm')
    }

    pvalue  <- getpvalue(u = cvm, eigen = ev)
    res     <- list(Statistic = cvm, pvalue = pvalue)

    return(res)

  } else if ( method == 'ad') {

    AD      <- getADStatistic(pit)

    if( gridpit ){
      ev    <- getEigenValues(S = Score, FI = fisher, pit, me = 'ad')
    }else{
      ev    <- getEigenValues_manualGrid(S = Score, FI = fisher, pit, M = ngrid, me = 'ad')
    }

    pvalue  <- getpvalue(u = AD, eigen = ev)
    res     <- list(Statistic = AD, pvalue = pvalue)

    return(res)

  }else{

    # Calculate both cvm and ad statistics

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
