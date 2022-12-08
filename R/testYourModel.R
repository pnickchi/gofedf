#' Apply Goodness of Fit (GoF) Test for user Specified Distribution Model
#'
#' @description This function can apply goodness of fit test based on empirical distribution function to any distribution defined by user.
#' In the case of no parameter estimation, the function requires score and probability transform inverse of data.
#' If there is any parameter that needs to be estimated, the function requires MLE of parameter too.
#'
#' @param x A numeric vector of length n for data points. The input must be numeric.
#' @param score can be a function provided by user that returns a matrix with n rows and p columns. The rows corresponds to each data point
#' and each column corresponds to each parameter in the model. If score is such a function, the first argument must take x (data points) and
#' the second argument must take a vector of mle estiamte of parameters. The function must return a matrix of n rows and p columns. testYourModel
#' function tries apply score function with mle as input to calculate score matrix.
#' The score input can instead be a matrix with n rows and p columns.
#' @param Fx can be a function to calculate the probability inverse transform of data.
#' A user provided function to calculate probability inverse transform of vector data, x. The first argument must be
#' vector of data values, x. The second argument must be a vector of unknown parameters. The number of parameters in this function
#' must match the number of parameters in score function (score_fn). See details for more examples.
#' @param mle Either null (if there is no parameter in the model) or a vector of mle estimate (s) of the unknown
#' parameter(s) in the model. The lenght of this vector must be the same as the number of columns of matrix
#' returned by score.
#' @param ngrid The number of equally spaced points to discritize the (0,1) interval to estimate the covariance of the stochastic process.
#' @param gridpit A Boolean indicator. If TRUE, ngrid is ignored and (0,1) interval is divided based on probability inverse transformed
#'  values. If FALSE (default value), (0,1) is divided into ngrid equally spaced points to estimate process.
#' @param precision The theory behind GoF based on empirical distribution function (edf) works well if the mle is indeed the root of
#' derivative of log likelihood. The user provided score function is evaluated at mle to verify its closeness to zero.
#' Precision is used to check how close score function is to zero. The default value of precision is 1e-6.
#' @param method a character to indicate which statistics to calculate.
#' The possible values are 'cvm' for Cramer-von Mises, 'ad' for Anderson-Darling, and 'both' for both methods. The default value is 'cvm'.
#' @return A list of two.
#' - Statistic: The Cramer-von-Mises statistic (and Anderson_Darling if method = 'both').
#' - pvalue: The approximate pvalue of the test. If method = 'both', two pvalues are returned.
#' @export
#'
#' @import stats
#'
#' @examples
#' set.seed(123)
#' # Example 1
#' # Generate some random data from Exponential dist
#' n <- 50
#' sim_data <- rexp(n, rate = 2)
#' # Estimate mle of scale parameter
#' theta    <- expMLE(x = sim_data)
#' testYourModel(x = sim_data, score = expScore, Fx = expFx, mle = theta, method = 'both')
#' # Example 2
#' # Generate some random data from Exponential dist
#' n <- 50
#' sim_data <- runif(n)
#' # Estimate mle of scale parameter
#' theta    <- expMLE(x = sim_data)
#' testYourModel(x = sim_data, score = expScore, Fx = expFx, mle = theta, method = 'cvm')
#' # Example 3
#' # Generate some random data from Normal dist
#' n <- 50
#' sim_data <- rnorm(n)
#' # Estimate mle of mean and sd
#' theta    <- normalMLE(x = sim_data)
#' testYourModel(x = sim_data, score = normalScore, Fx = normalPIT, mle = theta)
#' # Example 4
#' # Generate some random data from Normal dist
#' n <- 50
#' sim_data <- rnorm(n)
#' # Estimate mle of mean and sd
#' theta.value  <- normalMLE(x = sim_data)
#' score.matrix <- normalScore(x = sim_data, theta = theta.value)
#' pit.values   <- normalPIT(x = sim_data, theta = theta.value)
#' testYourModel(x = sim_data, score = score.matrix, Fx = pit.values, mle = theta.value)
testYourModel = function(x, score, Fx, mle = NULL, ngrid = length(x), gridpit = FALSE, precision = 1e-6, method = 'cvm'){

  if( !is.double(x) ){
    stop('x values must be numeric.')
  }

  if( length(x) < 2 ){
    stop('The numbers of data must be greater than or equal 2.')
  }

  if( !is.null(mle) ){

    if( any( c(!is.numeric(mle), !is.vector(mle)) ) ){
      stop('mle must be a numeric vector')
    }

  }

  check <- is.function(score) | is.matrix(score)
  if( !check ){
    stop('score must be either a function or a matrix.')
  }

  check <- is.function(Fx) | is.vector(Fx)
  if( !check ){
    stop('Fx must be either a function or a vector')
  }

  check <- is.function(score) & is.function(Fx)
  if( check ){

    # Get the number of parameters in score and Fx
    narg1 <- length( formals(score) ) - 1
    narg2 <- length( formals(Fx) )   - 1

    if( narg1 != narg2 ){
      stop('Number of parameters in score function and Fx function must be the same.')
    }

  }


  # Case 1: No parameter estimation
  if( is.null(mle) ){

    if( is.function(score) ){

      # Apply score function over vector data to get score function.
      score_matrix <- score(x)

      if( !is.matrix(score_matrix) ){
        stop('score function did not return a matrix.')
      }

    }else{
      score_matrix <- score
    }

    if( any(colSums(score_matrix) > precision) ){
      warning( paste0('Score matrix function is not zero at MLE. precision of ', precision, ' was used') )
    }

    if( is.function(Fx) ){
      # Apply Fx function to calculate probability inverse transform of data
      pit  <- lapply(x, FUN = Fx)
      pit  <- unlist(pit)
    }else{
      pit <- Fx
    }

    if( length(pit) != nrow(score_matrix) ){
      stop('number of rows in score matrix does not match the length of elements in pit vector.')
    }

    # Get the number of data points in sample
    n    <- length(x)

    # Calculate Fisher information matrix, estimate from sample score matrix
    fisher  <- (n-1)*var(score_matrix)/n

    # Get Eigen values
    if( gridpit ){
      ev    <- getEigenValues(S = score_matrix, FI = fisher, pit = pit, me = method)
    }else{
      ev    <- getEigenValues_manualGrid(S = score_matrix, FI = fisher, pit = pit, M = ngrid, me = method)
    }


    # Calculate Cramer-von-Mises, Anderson-Darling statistics or both

    if ( method == 'cvm') {

      # Calculate Cramer-von-Mises statistic
      U2        <- getCvMStatistic(pit)
      names(U2) <- 'Cramer-von-Mises Statistic'

      # Calculate pvalue
      pvalue  <- getpvalue(u = U2, eigen = ev)

      # Prepare a list to return statistic and pvalue
      res     <- list(Statistic = U2, pvalue = pvalue)

      return(res)

    } else if ( method == 'ad') {

      # Calculate Anderson-Darling statistic
      ad      <- getADStatistic(pit)
      names(ad) <- 'Anderson-Darling Statistic'

      # Calculate pvalue
      pvalue  <- getpvalue(u = ad, eigen = ev)

      # Prepare a list to return statistic and pvalue
      res     <- list(Statistic = ad, pvalue = pvalue)

      return(res)

    } else {

      # Calculate both cvm and ad statisitcs

      # 1. Do cvm calculation
      cvm        <- getCvMStatistic(pit)
      names(U2) <- 'Cramer-von-Mises Statistic'

      # Calculate pvalue
      cvm.pvalue  <- getpvalue(u = U2, eigen = ev)
      names(cvm.pvalue) <- 'pvalue for Cramer-von-Mises test'


      # 2. Do ad calculations
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

  # Case 2: Estimate parameters with MLE
  if( !is.null(mle) ){

    if( is.function(score) ){

      # Apply score function with mle over vector data to get score function.
      score_matrix <- score(x, mle)

      if( !is.matrix(score_matrix) ){
        stop('score function did not return a matrix.')
      }

    }else{
      score_matrix <- score
    }

    # Get the number of data points in sample
    n    <- length(x)

    # Get the number of parameters
    p <- length(mle)

    if( nrow(score_matrix) != n ){
      stop('The number of rows in score matrix do not match the sample size in x.')
    }

    if( ncol(score_matrix) != p){
      stop('The number of columns in score matrix do not match the number of estimated parameters in mle vector.')
    }

    if( any(colSums(score_matrix) > precision) ){
      warning( paste0('Score matrix function is not zero at MLE. precision of ', precision, ' was used') )
    }

    if( is.function(Fx) ){
      # Apply Fx function to calculate probability inverse transform of data
      pit  <- lapply(x, FUN = Fx, mle)
      pit  <- unlist(pit)
    }else{
      pit <- Fx
    }

    if( length(pit) != nrow(score_matrix) ){
      stop('number of rows in score matrix does not match the length of elements in pit vector.')
    }


    # Calculate Fisher information matrix, estimate from sample score matrix
    fisher  <- (n-1)*var(score_matrix)/n

    # Get Eigen values
    if( gridpit ){
      ev    <- getEigenValues(S = score_matrix, FI = fisher, pit = pit, me = method)
    }else{
      ev    <- getEigenValues_manualGrid(S = score_matrix, FI = fisher, pit = pit, M = ngrid, me = method)
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
      ad      <- getADStatistic(pit)
      names(ad) <- 'Anderson-Darling Statistic'

      # Calculate pvalue
      pvalue  <- getpvalue(u = ad, eigen = ev)

      # Prepare a list to return statistic and pvalue
      res     <- list(Statistic = ad, pvalue = pvalue)

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

}






