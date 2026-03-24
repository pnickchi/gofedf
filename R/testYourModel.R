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
#'   'both' to compute both cvm and ad statistics, and user for custom weight
#'   function. See weight_function for details about custom weight function.
#'
#' @param weight_function a function representing the weight function
#'   \eqn{w(u)} used to compute the weighted Cramér-von Mises statistic
#'   when \code{method = 'user'}. The function must take a numeric vector
#'   \eqn{u \in (0,1)} as input and return a numeric vector of the same
#'   length. The statistic is computed as
#'   \deqn{T_n = n \int_{0}^{1} w^2(u) \left( F_n(u) - u \right)^2 du}
#'   where \eqn{w^2(u)} is computed internally by squaring the supplied
#'   function. The default value is \code{NULL} and when \code{method} is
#'   \code{'cvm'}, \code{'ad'}, or \code{'both'} the weight_function is ignored.
#'   When \code{method = 'user'}, this argument must be provided, otherwise
#'   an error is returned.
#'
#' @return A list of two containing the following components:
#' - Statistic: the value of goodness-of-fit statistic.
#' - p-value: the approximate p-value for the goodness-of-fit test.
#'   if method = 'cvm' or method = 'ad', it returns a numeric value for the
#'   statistic and p-value. If method = 'both', it returns a numeric vector with
#'   two elements and one for each statistic. If method = 'user' it returns the
#'   weighted statistic.
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
#' scoreMatrix  <- IGScore(obs = sim_data, w = weights, mle = theta_hat)
#' pitvalues    <- IGPIT(obs = sim_data ,  w = weights, mle = theta_hat)
#'
#' # Apply the goodness-of-fit test.
#' testYourModel(pit = pitvalues, score = scoreMatrix)
#'
#' # Example to show custom weight function
#' w_cvm <- function(u) rep(1, length(u))
#' testYourModel(pit = pitvalues, score = scoreMatrix, method = 'user', weight_function = w_cvm)
#'
testYourModel = function(pit, score = NULL, discretize = FALSE, ngrid = length(pit), gridpit = TRUE, precision = 1e-9, method = 'cvm', weight_function = NULL){

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

  if( !(method %in% c('cvm','ad','both','user')) ){
    stop('method must be either cvm, ad, both, or user.')
  }

  if( method == 'user' & is.null(weight_function) ){
    stop('method is set to user but weight_function is NULL. You have to pass a weight function if you use user method.')
  }

  if( method == 'user' & !is.function(weight_function) ){
    stop('method is set to user but weight_function is not a valid function.')
  }

  if( !is.numeric(precision) ){
    stop('method must be numeric.')
  }

  if( method == 'user' & is.null(score) ){
    stop('Custom weight function in the case of no parameter estimation (score = NULL) is not supported in the current version of package.')
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

    }else if ( method == 'both'){

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

    }else{

      # Get the sample size
      n <- length(pit)

      # Define the identity matrix
      I <- diag( rep(1,n) )

      # Define matrix J
      J <- matrix(1/n, nrow = n, ncol = n)

      # Set the upper bound for integration - it is always constant
      UB <- n / (n + 1)

      # Create Q matrix
      Q <- matrix(0, n, n)

      # Loop over elements to fill Qij entry - use symmetry property
      for (i in 1:n){
        for (j in i:n){

          # Set the lower bound for this i and j
          LB <- max(i, j) / (n + 1)

          # Integrate
          Q[i, j] <- integrate(f = function(u) weight_function(u)^2, lower = LB, upper = UB, rel.tol = .Machine$double.eps^0.25)$value

          # Fill in the lower triangle
          Q[j, i] <- Q[i, j]
        }
      }

      # Define P matrix as P = (I - J) Q (I - J)
      P <- (I - J) %*% Q %*% (I - J)

      # Compute Eigen values
      ev      <- eigen(P, symmetric = TRUE, only.values = TRUE)$values

      # Compute weighted CvM statistics
      wcvm <- getWeightedStatistic(x = pit, w_function = weight_function)
      names(wcvm) <- 'Weighted Cramer-von-Mises Statistic'

      # Calculate pvalue
      pvalue  <- getpvalue(u = wcvm, eigen = ev)

      # Prepare a list to return statistic and pvalue
      res     <- list(Statistic = wcvm, pvalue = pvalue)

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
      warning( paste0('score is not zero at the MLE that you used. \n precision of ', precision, ' was used. \n The asymptotic assumptions may not be accurate.') )
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

      if( method == 'user' ){

        # Compute P matrix
        P <- computeMatrix(n, score, method = method, w_function = weight_function)

        # Adjust for the number of estimated parameters
        P <- P / (n-ncol(score)-1)

        # Compute eigenvalues
        ev <- eigen(P, only.values = TRUE, symmetric = TRUE)$values

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

      }else if ( method == 'both' ){

        # Computue cvm and AD statistics along with their pvalues
        cvm <- getCvMStatistic(pit)
        cvm.pvalue  <- getpvalue(u = cvm, eigen = ev_cvm)

        AD  <- getADStatistic(pit)
        ad.pvalue  <- getpvalue(u = AD, eigen = ev_ad)

        gof.stat        <- c(cvm, AD)
        names(gof.stat) <- c('Cramer-von-Mises Statistic','Anderson-Darling Statistic')

        # Prepare a list to return statistic and pvalue
        res     <- list(Statistics = gof.stat, pvalue = c(cvm.pvalue, ad.pvalue) )

        return(res)

      }else{

        # Compute weighted CvM statistics
        wcvm <- getWeightedStatistic(x = pit, w_function = weight_function)
        names(wcvm) <- 'Weighted Cramer-von-Mises Statistic'

        # Calculate pvalue
        pvalue  <- getpvalue(u = wcvm, eigen = ev)

        # Prepare a list to return statistic and pvalue
        res     <- list(Statistic = wcvm, pvalue = pvalue)

        return(res)

      }

    }

    #
    # Use the estimated covariance function and turning integral equation into a matrix equation.
    # This is where discrritize = TRUE is applied
    #

    # Calculate Fisher information matrix by computing the variance of score from the sample.
    fisher  <- (n-1)*var(score)/n


    # Calculate Cramer-von-Mises, Anderson-Darling statistics or both
    if ( method == 'cvm' ) {

      # Calculate Cramer-von-Mises statistic
      cvm        <- getCvMStatistic(pit)
      names(cvm) <- 'Cramer-von-Mises Statistic'

      # Get Eigen values
      if( gridpit ){
        ev    <- getEigenValues(S = score, FI = fisher, pit = pit, me = method)
      }else{
        ev    <- getEigenValues_manualGrid(S = score, FI = fisher, pit = pit, M = ngrid, me = method)
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
        ev    <- getEigenValues(S = score, FI = fisher, pit = pit, me = method)
      }else{
        ev    <- getEigenValues_manualGrid(S = score, FI = fisher, pit = pit, M = ngrid, me = method)
      }

      # Calculate pvalue
      pvalue  <- getpvalue(u = AD, eigen = ev)

      # Prepare a list to return statistic and pvalue
      res     <- list(Statistic = AD, pvalue = pvalue)

      return(res)

    }else if ( method == 'both' ){

      # Calculate for both cvm and ad statistics

      # 1. Do cvm calculation

      # Compute Eigen values for cvm
      if( gridpit ){
        ev    <- getEigenValues(S = score, FI = fisher, pit, me = 'cvm')
      }else{
        ev    <- getEigenValues_manualGrid(S = score, FI = fisher, pit, M = ngrid, me = 'cvm')
      }

      cvm        <- getCvMStatistic(pit)
      names(cvm) <- 'Cramer-von-Mises Statistic'

      # Calculate pvalue
      cvm.pvalue  <- getpvalue(u = cvm, eigen = ev)
      names(cvm.pvalue) <- 'pvalue for Cramer-von-Mises test'


      # 2. Do ad calculations

      # Compute Eigen values for ad
      if( gridpit ){
        ev    <- getEigenValues(S = score, FI = fisher, pit, me = 'ad')
      }else{
        ev    <- getEigenValues_manualGrid(S = score, FI = fisher, pit, M = ngrid, me = 'ad')
      }

      ad      <- getADStatistic(pit)
      names(ad) <- 'Anderson-Darling Statistic'

      # Calculate pvalue
      ad.pvalue  <- getpvalue(u = ad, eigen = ev)
      names(ad.pvalue) <- 'Anderson-Darling test'

      # Prepare a list to return both statistics and their approximate pvalue
      res     <- list(Statistics = c(cvm, ad), pvalue = c(cvm.pvalue, ad.pvalue) )
      return(res)

    }else{

      # Compute Eigen values
      if( gridpit ){
        ev    <- getEigenValues(S = score, FI = fisher, pit, me = method, w_function = weight_function)
      }else{
        ev    <- getEigenValues_manualGrid(S = score, FI = fisher, pit, M = ngrid, me = method, w_function = weight_function)
      }

      # Compute Cramer-von-Mises statistic
      wcvm      <- getWeightedStatistic(x = pit, w_function = weight_function)
      names(wcvm) <- 'Weighted Cramer-von-Mises Statistic'

      # Compute pvalue
      pvalue  <- getpvalue(u = wcvm, eigen = ev)
      res     <- list(Statistic = wcvm, pvalue = pvalue)

      return(res)

    }

  }

}
