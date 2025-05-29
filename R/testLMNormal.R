#' Apply Goodness of Fit Test to Residuals of a Linear Model
#'
#' @description \code{testLMNormal} is used to check the normality assumption of
#'   residuals in a linear model. This function can take the response variable
#'   and design matrix, fit a linear model, and apply the goodness-of-fit test.
#'   Conveniently, it can take an object of class "lm" and directly applies the
#'   goodness-of-fit test. The function returns a goodness-of-fit statistic
#'   along with an approximate p-value.
#'
#' @param x is either a numeric vector or a design matrix. In the design matrix,
#'   rows indicate observations and columns presents covariates.
#'
#' @param y is a vector of numeric values with the same number of observations
#'   or number of rows as x.
#'
#' @param fit an object of class "lm" returned by  \code{\link{lm}} function in
#'   \code{\link{stats}} package. The default value of fit is NULL. If any
#'   object is provided, \code{x} and \code{y} will be ignored and the class of
#'   object is checked. If you pass an object to \code{fit} make sure to return
#'   the design matrix by setting \code{x} = \code{TRUE} and the response
#'   variable by setting in \code{y} = \code{TRUE} in \code{\link{lm}} function.
#'   To read more about this see the help documentation for \code{\link{lm}}
#'   function or see the example below.
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
#' set.seed(123)
#' n <- 50
#' p <- 5
#' x <- matrix( runif(n*p), nrow = n, ncol = p)
#' e <- rnorm(n)
#' b <- runif(p)
#' y <- x %*% b + e
#' testLMNormal(x, y)
#' # Or pass lm.fit object directly:
#' lm.fit <- lm(y ~ x, x = TRUE, y = TRUE)
#' testLMNormal(fit = lm.fit)
testLMNormal = function(x, y, fit = NULL, discretize = FALSE, ngrid = length(y), gridpit = TRUE, hessian = FALSE, method = 'cvm'){


  if( is.null(fit) ){


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
      stop('Method must be either cvm, ad, or both.')
    }

    if( is.vector(x) ){
      n           <- length(x)
      int         <- rep(1, n)
      x           <- cbind(int,x)
      colnames(x) <- c('Intercept', 'x')
    }

    if( is.matrix(x) ){
      n    <- nrow(x)
      int  <- rep(1, n)
      x    <- cbind(int,x)
    }

  }

  if( !is.null(fit) ){

    if (!inherits(fit, 'lm')){
      stop('The fit must be \'lm\' object.')
    }

    if( !is.matrix(fit$x) ){
      stop('fit object must have the design matrix corresponding to the model. \n Consider setting x = TRUE in lm function to return x matrix.')
    }

    if( !is.vector(fit$y) ){
      stop('fit object must contain the response variable. \n Consider setting y = TRUE in lm function to return reponse variable.')
    }

    x <- fit$x
    y <- fit$y
    n <- length(y)

  }

  # Apply linear model with normal assumption
  temp   <- applyLMNormal(x = x, y = y)

  # Extract score function, pit values, and MLE
  Score  <- temp$Score
  pit    <- temp$pit
  par    <- temp$par

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
    fisher <- lmFisherByHessian(x = x, y = y, theta = par)
  }else{
    fisher <- (n-1)*var(Score)/n
  }

  # Compute Eigen values
  if( gridpit ){
    ev    <- getEigenValues(S = Score, FI = fisher, pit, me = method)
  }else{
    ev    <- getEigenValues_manualGrid(S = Score, FI = fisher, pit, M = ngrid, me = method)
  }

  if( method == 'cvm'){

    # Compute Cramer-von-Mises statistic
    cvm      <- getCvMStatistic(pit)

    # Compute pvalue
    pvalue  <- getpvalue(u = cvm, eigen = ev)
    res     <- list(Statistic = cvm, pvalue = pvalue)

    return(res)

  } else if ( method == 'ad') {

    # Compute Anderson-Darling statistic
    AD      <- getADStatistic(pit)

    # Compute pvalue
    pvalue  <- getpvalue(u = AD, eigen = ev)

    res     <- list(Statistic = AD, pvalue = pvalue)
    return(res)

  }else{

    # Calculate both cvm and ad statistics

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
