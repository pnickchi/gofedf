#' Apply Goodness of Fit (GoF) Test for a Specified Distribution Model
#'
#' @description This function can apply goodness of fit test based on empirical distribution function to any distribution defined by user.
#' The function requires score and probability transform inverse of data.
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
#' @param mle A vector of mle estimate (s) of the unknown parameter(s) in the model. The lenght of this vector must be the same as the number of
#' columns of matrix returned by score.
#' @param ngrid The number of equally spaced points to discritize the (0,1) interval to estimate the covariance of the stochastic process.
#' @param gridpit A Boolean indicator. If TRUE, ngrid is ignored and (0,1) interval is divided based on probability inverse transformed
#'  values. If FALSE (default value), (0,1) is divided into ngrid equally spaced points to estimate process.
#' @param precision The theory behind GoF based on empirical distribution function (edf) works well if the mle is indeed the root of
#' derivative of log likelihood. The user provided score function is evaluated at mle to verify its closeness to zero.
#' Precision is used to check how close score function is to zero. The default value of precision is 1e-6.
#' @param method a character string to indicate which statistics to calculate.
#' The possible values are 'cvm' for Cramer-von Mises and 'ad' for Anderson-Darling.
#' The default value is 'cvm'.
#' @return A list of two.
#' - Statistic: The Cramer-von-Mises statistic.
#' - pvalue: The approximated pvalue for the GoF test based on EDF.
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
#' testYourModel(x = sim_data, score = expScore, Fx = expFx, mle = theta)
#' # Example 2
#' # Generate some random data from Exponential dist
#' n <- 50
#' sim_data <- runif(n)
#' # Estimate mle of scale parameter
#' theta    <- expMLE(x = sim_data)
#' testYourModel(x = sim_data, score = expScore, Fx = expFx, mle = theta)
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
testYourModel = function(x, score, Fx, mle, ngrid = length(x), gridpit = FALSE, precision = 1e-6, method = 'cvm'){

  if( !is.double(x) ){
    stop('x values must be numeric.')
  }

  if( length(x) < 2 ){
    stop('The numbers of data must be greater than or equal 2.')
  }

  if( any( c(!is.numeric(mle), !is.vector(mle)) ) ){
    stop('mle must be a numeric vector')
  }

  check <- is.function(score) | is.matrix(score)
  if( !check ){
    stop('score must be either a function or a matrix.')
  }

  check <- is.function(Fx) | is.vector(Fx)
  if( !check ){
    stop('Fx must be either a function or a vector')
  }

  if( is.function(score) & is.function(Fx) ){

    # Get the number of data points in sample
    n    <- length(x)

    # Get the number of parameters
    p <- length(mle)

    # Get the number of parameters in score and Fx
    narg1 <- length( formals(score) ) - 1
    narg2 <- length( formals(Fx) )   - 1

    if( narg1 != narg2 ){
      stop('Number of parameters in score function and Fx function must be the same')
    }

    # Apply score function over vector data to get score function.
    score_matrix <- score(x, mle)

    if( !is.matrix(score_matrix) ){
      stop('score function did not return a matrix.')
    }

    if( any(colSums(score_matrix) > precision) ){
      warning( paste0('Score matrix function is not zero at MLE. precision of ', precision, ' was used') )
    }

    # Apply Fx function to calculate probability inverse transform of data
    pit  <- lapply(x, FUN = Fx, mle)
    pit  <- unlist(pit)

    fisher  <- (n-1)*var(score_matrix)/n

    if( gridpit ){
      ev    <- getEigenValues(S = score_matrix, FI = fisher, pit = pit, me = method)
    }else{
      ev    <- getEigenValues_manualGrid(S = score_matrix, FI = fisher, pit = pit, M = ngrid, me = method)
    }


    U2      <- getCvMStatistic(pit)
    pvalue  <- getpvalue(u = U2, eigen = ev)

    res     <- list(Statistic = U2, pvalue = pvalue)
    return(res)

  }

  if( is.matrix(score)   & is.vector(Fx) ){

    # Get the number of data points in sample
    n    <- length(x)

    # Get the number of parameters
    p <- length(mle)

    if( nrow(score) != n ){
      stop('The number of rows in score matrix do not match the sample size in x.')
    }

    if( ncol(score) != p){
      stop('The number of columns in score matrix do not match the number of estimated parameters in mle vector.')
    }

    if( any( is.na(x) ) ){
      stop('There are missing values in vector x.')
    }

    if( any( is.na(score) ) ){
      stop('There are missing values in matrix score')
    }

    if( any( is.na(Fx) ) ){
      stop('There are missing values in Fx values.')
    }

    if( any(colSums(score) > precision) ){
      warning( paste0('Score matrix function is not zero at MLE. precision of ', precision, ' was used') )
    }

    fisher  <- (n-1)*var(score)/n

    if( gridpit ){
      ev    <- getEigenValues(S = score, FI = fisher, pit = Fx, me = method)
    }else{
      ev    <- getEigenValues_manualGrid(S = score, FI = fisher, pit = Fx, M = ngrid, me = method)
    }

    U2      <- getCvMStatistic(x = Fx)
    pvalue  <- getpvalue(u = U2, eigen = ev)

    res     <- list(Statistic = U2, pvalue = pvalue)
    return(res)

  }

}
