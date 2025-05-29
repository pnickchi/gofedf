#' Apply Goodness of Fit Test for Normal Distribution
#'
#' @description Performs the goodness-of-fit test based on empirical
#'   distribution function to check if an i.i.d sample follows a Normal
#'   distribution.
#'
#' @param x a non-empty numeric vector of sample data.
#'
#' @param discretize If \code{TRUE}, the covariance function of \eqn{W_{n}(u)}
#'   process is evaluated at some data points (see \code{ngrid} and
#'   \code{gridpit}), and the integral equation is replaced by a matrix
#'   equation. If \code{FALSE} (the default value), the covariance function is
#'   first estimated, and then the integral equation is solved to find the
#'   eigenvalues. The results of our simulations recommend using the estimated
#'   covariance for solving the integral equation. The parameters \code{ngrid},
#'   \code{gridpit}, and \code{hessian} are only relevant when \code{discretize
#'   = TRUE}.
#'
#' @param ngrid the number of equally spaced points to discretize the (0,1)
#'   interval for computing the covariance function.
#'
#' @param gridpit logical. If \code{TRUE} (the default value), the parameter
#'   ngrid is ignored and (0,1) interval is divided based on probability
#'   integral transforms or PITs obtained from the sample. If \code{FALSE}, the
#'   interval is divided into ngrid equally spaced points for computing the
#'   covariance function.
#'
#' @param hessian logical. If \code{TRUE} the Fisher information matrix is
#'   estimated by the observed Hessian Matrix based on the sample. If
#'   \code{FALSE} (the default value) the Fisher information matrix is estimated
#'   by the variance of the observed score matrix.
#'
#' @param method a character string indicating which goodness-of-fit statistic
#'   is to be computed. The default value is 'cvm' for the Cramer-von-Mises
#'   statistic. Other options include 'ad' for the Anderson-Darling statistic,
#'   and 'both' to compute both cvm and ad.
#'
#' @return A list of two containing the following components:
#' - Statistic: the value of goodness-of-fit statistic.
#' - p-value: the approximate p-value for the goodness-of-fit test.
#'   if method = 'cvm' or method = 'ad', it returns a numeric value for the
#'   statistic and p-value. If method = 'both', it returns a numeric vector with
#'   two elements and one for each statistic.
#'
#' @export
#'
#' @examples
#'
#' set.seed(123)
#' sim_data <- rnorm(n = 50)
#' testNormal(x = sim_data)
#' sim_data <- rgamma(n = 50, shape = 3)
#' testNormal(x = sim_data)
testNormal = function(x, discretize = FALSE, ngrid = length(x), gridpit = TRUE, hessian = FALSE, method = 'cvm'){

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

  if( anyNA(x) ){
    x <- x[ !is.na(x) ]
    warning('NA found in x and automatically removed.')
  }

  # Get the sample size
  n       <- length(x)

  # Apply Normal distribution
  temp    <- applyNormal(x)

  # Compute score function, pit values, and MLE
  Score   <- temp$Score
  pit     <- temp$pit
  par     <- temp$par


  # Use the estimated covariance function when solving the integral equation
  if(!discretize){

    # Find the rank of sorted pits
    sort_indx <- order(pit)

    # Reorder the rows of score matrix according to the ranks in pit
    Score     <- Score[sort_indx,]

    if( method == 'cvm' | method == 'ad' ){

      # Compute P matrix
      P <- computeMatrix(n, Score, method = method)

      # Adjust for the number of estimated parameters
      P <- P / (n-ncol(Score)-1)

      # Compute eigenvalues
      ev <- eigen(P, only.values = TRUE, symmetric = TRUE)$values

    }

    if( method == 'both' ){

      # Compute P matrix, adjust for number of estimated parameters, and compute eigenvalues for the case of cvm
      P_cvm  <- computeMatrix(n, Score, method = 'cvm')
      P_cvm  <- P_cvm/(n-ncol(Score)-1)
      ev_cvm <- eigen(P_cvm, only.values = TRUE, symmetric = TRUE)$values

      # Compute P matrix, adjust for number of estimated parameters, and compute eigenvalues for the case of ad
      P_ad  <- computeMatrix(n, Score, method = 'ad')
      P_ad  <- P_ad/(n-ncol(Score)-1)
      ev_ad <- eigen(P_ad, only.values = TRUE, symmetric = TRUE)$values

    }


    # Compute gof statistics and pvalue according to the requested method
    if( method == 'cvm' ){

      # Compute CvM statistics
      cvm <- getCvMStatistic(pit)
      names(cvm) <- 'Cramer-von-Mises Statistic'

      # Calculate pvalue
      pvalue  <- getpvalue(u = cvm, eigen = ev)

      # Prepare a list to return statistic and pvalue
      res     <- list(Statistic = cvm, pvalue = pvalue)

      return(res)

    } else if ( method == 'ad' ){

      AD <- getADStatistic(pit)
      names(AD) <- 'Anderson-Darling Statistic'

      # Calculate pvalue
      pvalue  <- getpvalue(u = AD, eigen = ev)

      # Prepare a list to return statistic and pvalue
      res     <- list(Statistic = AD, pvalue = pvalue)

      return(res)

    }else{

      # Fix this oart needs to have two diff sets of ev
      cvm <- getCvMStatistic(pit)
      cvm.pvalue  <- getpvalue(u = cvm, eigen = ev_cvm)

      AD  <- getADStatistic(pit)
      ad.pvalue  <- getpvalue(u = AD, eigen = ev_ad)

      gof.stat        <- c(cvm, AD)
      names(gof.stat) <- c('Cramer-von-Mises Statistic','Anderson-Darling Statistic')

      # Prepare a list to return statistic and pvalue
      res     <- list(Statistics = gof.stat, pvalue = c(cvm.pvalue, ad.pvalue) )

      return(res)

    }

  }



  # Lines below here are used when there is discritization to compute the covariance of W_{n}(u) process.


  # Compute Fisher information matrix
  if( hessian ){
    fisher <- normalFisherByHessian(par)
  }else{
    fisher <- (n-1)*var(Score)/n
  }



  if( method == 'cvm'){

    # Compute Eigen values
    if( gridpit ){
      ev    <- getEigenValues(S = Score, FI = fisher, pit, me = method)
    }else{
      ev    <- getEigenValues_manualGrid(S = Score, FI = fisher, pit, M = ngrid, me = method)
    }

  # Compute Cramer-von-Mises statistic
    cvm      <- getCvMStatistic(pit)

    # Compute pvalue
    pvalue  <- getpvalue(u = cvm, eigen = ev)
    res     <- list(Statistic = cvm, pvalue = pvalue)

    return(res)

  } else if ( method == 'ad') {

    # Compute Eigen values
    if( gridpit ){
      ev    <- getEigenValues(S = Score, FI = fisher, pit, me = method)
    }else{
      ev    <- getEigenValues_manualGrid(S = Score, FI = fisher, pit, M = ngrid, me = method)
    }

    # Compute Anderson-Darling statistic
    AD      <- getADStatistic(pit)

    # Compute pvalue
    pvalue  <- getpvalue(u = AD, eigen = ev)

    res     <- list(Statistic = AD, pvalue = pvalue)
    return(res)

  }else{

    # Calculate for both cvm and ad statistics

    # 1. Do cvm calculation

    # Compute Eigen values for cvm
    if( gridpit ){
      ev    <- getEigenValues(S = Score, FI = fisher, pit, me = 'cvm')
    }else{
      ev    <- getEigenValues_manualGrid(S = Score, FI = fisher, pit, M = ngrid, me = 'cvm')
    }

    cvm        <- getCvMStatistic(pit)
    names(cvm) <- 'Cramer-von-Mises Statistic'

    # Calculate pvalue
    cvm.pvalue  <- getpvalue(u = cvm, eigen = ev)
    names(cvm.pvalue) <- 'pvalue for Cramer-von-Mises test'


    # 2. Do ad calculations

    # Compute Eigen values for ad
    if( gridpit ){
      ev    <- getEigenValues(S = Score, FI = fisher, pit, me = 'ad')
    }else{
      ev    <- getEigenValues_manualGrid(S = Score, FI = fisher, pit, M = ngrid, me = 'ad')
    }

    ad      <- getADStatistic(pit)
    names(ad) <- 'Anderson-Darling Statistic'

    # Calculate pvalue
    ad.pvalue  <- getpvalue(u = ad, eigen = ev)
    names(ad.pvalue) <- 'Anderson-Darling test'

    # Prepare a list to return both statistics and their approximate pvalue
    res     <- list(Statistics = c(cvm, ad), pvalue = c(cvm.pvalue, ad.pvalue) )
    return(res)

  }

}
