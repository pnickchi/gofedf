#' Apply the Goodness of Fit Test Based on Empirical Distribution Function to
#' Any Likelihood Model.
#'
#' @description This function applies the goodness-of-fit test based on
#'   empirical distribution function. It requires certain inputs depending on
#'   whether the model involves parameter estimation or not. If the model is
#'   known and there is no parameter estimation, the function requires the
#'   probability transformed (or pit) values of the sample. This ought to be a
#'   numeric vector. If there is parameter estimation in the model, the function
#'   additionally requires the score as a matrix with n rows and p columns,
#'   where n is the sample size and p is the number of estimated parameters. The
#'   function checks if the sum of columns in score is near zero at the
#'   estimated parameter (which is assumed to be the maximum likelihood
#'   estimate).
#'
#' @param pit The probability transformed (or pit) values of the sample which
#'   ought to be a numeric vector.
#'
#' @param score The default value is null and refers to no parameter estimation
#'   case. If there is parameter estimation, the score must be a matrix with n
#'   rows and p columns, where n is the sample size and p is the number of
#'   estimated parameters.
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
#' @param ngrid The number of equally spaced points to discretize the
#'   (0,1)interval for computing the covariance function.
#'
#' @param gridpit logical. If \code{TRUE} (the default value), the parameter
#'   ngrid is ignored and (0,1) interval is divided based on probability
#'   integral transforms or PITs obtained from the sample. If \code{FALSE}, the
#'   interval is divided into ngrid equally spaced points for computing the
#'   covariance function.
#'
#' @param precision The theory behind goodness-of-fit test based on empirical
#'   distribution function (edf) works well if the MLE is indeed the root of
#'   derivative of log likelihood function. A precision of 1e-9 (default value)
#'   is used to check this. A warning message is generated if the score
#'   evaluated at MLE is not close enough to zero.
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
#' theta_hat    <- IGMLE(obs = sim_data,   w = weights)
#' ScoreMatrix  <- IGScore(obs = sim_data, w = weights, mle = theta_hat)
#' pitvalues    <- IGPIT(obs = sim_data ,  w = weights, mle = theta_hat)
#'
#' # Apply the goodness-of-fit test.
#' testYourModel(pit = pitvalues, score = ScoreMatrix)
#'
testYourModel = function(pit, score = NULL, discretize = FALSE, ngrid = length(pit), gridpit = TRUE, precision = 1e-9, method = 'cvm'){

  if( !is.numeric(pit) ){
    stop('pit values must be numeric.')
  }

  if( !is.vector(pit) ){
    stop('pit must be a vector of numeric values.')
  }

  if( length(pit) < 2 ){
    stop('The number of observations in pit must be greater than or equal two.')
  }

  if( any(pit < 0) | any(pit > 1) ){
    stop('pit values must be between zero and one.')
  }

  if( !is.logical(gridpit) ){
    stop('gridpit must be either TRUE or FALSE.')
  }

  if( !is.logical(discretize) ){
    stop('discretize must be either TRUE or FALSE.')
  }

  if( !is.vector(method) | length(method) > 1){
    stop('method must be a character string with length one.')
  }

  if( !(method %in% c('cvm','ad','both')) ){
    stop('method must be either cvm, ad, or both.')
  }

  if( !is.numeric(precision) ){
    stop('method must be numeric.')
  }

  if( anyNA(pit) ){
    pit <- pit[ !is.na(pit) ]

    if( !is.null(score) ){
       indx  <- which(is.na(pit))
       score <- score[-indx,]
    }

    warning('NA found in pit and automatically removed. \n Corresponding rows in score were also removed.')
  }


  #
  # Case 1: if score is null then it means there was no parameter estimation in the model.
  # The model is fully specified and the covariance function of the  stochastic process,
  # $\hat W_{n}(u)$, is simply min(s,t)-st for 0 <= s,t <=1. The Eigen values are:
  # (i) for Cramer-von-Mises:  1 / ( pi^2 * j^2 ); j=1,2,3...
  # (ii) for Anderson-Darling: 1 / ( j * (j+1) ); j=1,2,3,...
  #
  if( is.null(score) ){

    # Compute the 100 leading Eigenvalues of the covariance function
    j  <- 1:100

    if ( method == 'cvm') {

      ev <- 1 / ( pi^2 * j^2 )

      # Calculate Cramer-von-Mises statistic
      cvm        <- getCvMStatistic(pit)
      names(cvm) <- 'Cramer-von-Mises Statistic'

      # Calculate pvalue
      pvalue  <- getpvalue(u = cvm, eigen = ev)

      # Prepare a list to return statistic and pvalue
      res     <- list(Statistic = cvm, pvalue = pvalue)

      return(res)

    } else if ( method == 'ad') {

      ev <- 1 / ( j * (j+1) )

      # Calculate Anderson-Darling statistic
      AD        <- getADStatistic(pit)
      names(AD) <- 'Anderson-Darling Statistic'

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

      # Calculate pvalue
      ev <- 1 / ( pi^2 * j^2 )
      cvm.pvalue  <- getpvalue(u = cvm, eigen = ev)
      names(cvm.pvalue) <- 'pvalue for Cramer-von-Mises test'

      # 2. Do AD calculations
      ad      <- getADStatistic(pit)
      names(ad) <- 'Anderson-Darling Statistic'

      # Calculate pvalue
      ev <- 1 / ( j * (j+1) )
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
  # The covariance function of the stochastic process, $\hat W_{n}(u)$, needs to be calculated.
  # We estimate this covariance function.
  #
  if( !is.null(score) ){

    # Check if the score is a matrix
    if( !is.matrix(score) ){
      stop('score must be a matrix.')
    }

    # Get the length of pit.
    n    <- length(pit)

    # Check if score dimension is correct
    if( nrow(score) != n ){
      stop('The number of rows in score matrix do not match the sample size in pit.')
    }

    # Check if score is zero at MLE
    if( any(colSums(score) > precision) ){
      warning( paste0('Score is not zero at the MLE that you used. \n precision of ', precision, ' was used. \n The asymptotic assumptions may not be accurate.') )
    }

    #
    # Use the estimated covariance function and use the analytical solution of integral equation to find eigen values.
    #
    if(!discretize){

      # Find the rank of sorted pits
      sort_indx <- order(pit)

      # Reorder the rows of score matrix according to the ranks in pit
      score     <- score[sort_indx,]

      if( method == 'cvm' | method == 'ad' ){

        # Compute P matrix
        P <- computeMatrix(n, score, method = method)

        # Adjust for the number of estimated parameters
        P <- P / (n-ncol(score)-1)

        # Compute eigenvalues
        ev <- eigen(P, only.values = TRUE, symmetric = TRUE)$values

      }

      if( method == 'both' ){

        # Compute P matrix, adjust for number of estimated parameters, and compute eigenvalues for the case of cvm
        P_cvm  <- computeMatrix(n, score, method = 'cvm')
        P_cvm  <- P_cvm/(n-ncol(score)-1)
        ev_cvm <- eigen(P_cvm, only.values = TRUE, symmetric = TRUE)$values

        # Compute P matrix, adjust for number of estimated parameters, and compute eigenvalues for the case of ad
        P_ad  <- computeMatrix(n, score, method = 'ad')
        P_ad  <- P_ad/(n-ncol(score)-1)
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


    #
    # Use the estimated covariance function and turning integral equation into a matrix equation.
    # This is where discrritize = TRUE is applied
    #

    # Calculate Fisher information matrix by computing the variance of score from the sample.
    fisher  <- (n-1)*var(score)/n

    # Get Eigen values
    if( gridpit ){
      ev    <- getEigenValues(S = score, FI = fisher, pit = pit, me = method)
    }else{
      ev    <- getEigenValues_manualGrid(S = score, FI = fisher, pit = pit, M = ngrid, me = method)
    }


    # Calculate Cramer-von-Mises, Anderson-Darling statistics or both
    if ( method == 'cvm' ) {

      # Calculate Cramer-von-Mises statistic
      cvm        <- getCvMStatistic(pit)
      names(cvm) <- 'Cramer-von-Mises Statistic'

      # Calculate pvalue
      pvalue  <- getpvalue(u = cvm, eigen = ev)

      # Prepare a list to return statistic and pvalue
      res     <- list(Statistic = cvm, pvalue = pvalue)

      return(res)

    } else if ( method == 'ad') {

      # Calculate Anderson-Darling statistic
      AD      <- getADStatistic(pit)
      names(AD) <- 'Anderson-Darling Statistic'

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

      # Calculate pvalue
      cvm.pvalue  <- getpvalue(u = cvm, eigen = ev)
      names(cvm.pvalue) <- 'pvalue for Cramer-von-Mises test'


      # 2. Do ad calculations
      AD        <- getADStatistic(pit)
      names(AD) <- 'Anderson-Darling Statistic'

      # Calculate pvalue
      ad.pvalue  <- getpvalue(u = AD, eigen = ev)
      names(ad.pvalue) <- 'Anderson-Darling test'


      # Prepare a list to return both statistics and their approximate pvalue
      res     <- list(Statistics = c(cvm, AD), pvalue = c(cvm.pvalue, ad.pvalue) )
      return(res)

    }

  }

}
