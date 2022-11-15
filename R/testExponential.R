#' Apply Goodness of Fit Test for Exponential Distribution
#'
#' @param x A numeric vector of data values with length n.
#' @param ngrid The number of equally spaced points to discritize the (0,1) interval to estimate the covariance of the process.
#' @param gridpit A Boolean indicator. If TRUE, ngrid is ignored and (0,1) interval is divided based on probability inverse transformed
#'  values. If FALSE (default value), (0,1) is divided into ngrid equally spaced points to estimate process.
#' @param hessian A Boolean indicator. If TRUE, the Fisher information matrix is estimated by observed Hessian Matrix.
#' If False (the default value), the Fisher information matrix is estimated by observed score function.
#'
#' @param method a character string to indicate which statistics to calculate.
#' The possible values are 'cvm' for Cramer-von Mises and 'ad' for Anderson-Darling.
#' The default value is 'cvm'.

#' @return A list of two.
#' - Statistic: The Cramer-von-Mises statistic.
#' - pvalue: The approximated pvalue for the GoF test based on EDF.
#' @export
#'
#' @examples
#' set.seed(123)
#' sim_data <- rexp(n = 100)
#' testExponential(x = sim_data)
testExponential = function(x, ngrid = length(x), gridpit = FALSE, hessian = FALSE, method = 'cvm'){

  if( any( x < 0 ) ){
    stop('x values must be positive for Exponential.')
  }

  n       <- length(x)
  temp    <- applyExponential(x)
  Score   <- temp$Score
  pit     <- temp$pit
  par     <- temp$par

  if( hessian ){
    fisher <- matrix( n / (par[1])^2, nrow = 1, ncol = 1 )
  }else{
    fisher <- (n-1)*var(Score)/n
  }

  if( gridpit ){
    ev    <- getEigenValues(S = Score, FI = fisher, pit, me = method)
  }else{
    ev    <- getEigenValues_manualGrid(S = Score, FI = fisher, pit, M = ngrid, me = method)
  }

  U2      <- getCvMStatistic(pit)
  pvalue  <- getpvalue(u = U2, eigen = ev)

  res     <- list(Statistic = U2, pvalue = pvalue)
  return(res)

}
