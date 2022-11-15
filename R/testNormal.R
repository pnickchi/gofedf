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

  if( !is.integer(ngrid) | !(ngrid > 0) ){
    stop('ngrid must be a positive number.')
  }

  if( !is.logical(gridpit) ){
    stop('gridpit must be a logical value.')
  }

  if( !is.logical(hessian) ){
    stop('hessian must be a logical value.')
  }

  if( !is.vector(method) | length(method) > 1){
    stop('method must be a character string.')
  }

  if( !(method %in% c('cvm','ad')) ){
     stop('method must be either cvm or ad.')
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

  if( gridpit ){
    ev    <- getEigenValues(S = Score, FI = fisher, pit, me = method)
  }else{
    ev    <- getEigenValues_manualGrid(S = Score, FI = fisher, pit, M = ngrid, me = method)
  }

  if( method == 'cvm'){
    U2      <- getCvMStatistic(pit)
    pvalue  <- getpvalue(u = U2, eigen = ev)
    res     <- list(Statistic = U2, pvalue = pvalue)
  }else{
    AD      <- getADStatistic(pit)
    pvalue  <- getpvalue(u = AD, eigen = ev)
    res     <- list(Statistic = AD, pvalue = pvalue)
  }

  return(res)

}
