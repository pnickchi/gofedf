#' Apply Goodness of Fit Test to Residuals of a Linear Model
#'
#' @description \code{testLMNormal} is used to check the normality assumption of residuals in a linear model. This function
#' can take the response variable and design matrix, fit a linear model, and apply the goodness-of-fit test. Conveniently,
#' it can take an object of class "lm" and directly applies the goodness-of-fit test. The function returns a goodness-of-fit
#' statistic along with an approximate pvalue.
#'
#' @param x is either a numeric vector or a design matrix. In the design matrix, rows indicate observations and columns
#' presents covariats.
#'
#' @param y is a vector of numeric values with the same number of observations or number of rows as x.
#'
#' @param fit an object of class "lm" returned by  \code{\link{lm}} function in \code{\link{stats}} package. The default value of
#' fit is NULL. If any object is provided, \code{x} and \code{y} will be ignored and the class of object is checked. If you pass
#' an object to \code{fit} make sure to return the design matrix by setting \code{x} = \code{TRUE} and the response variable by
#' setting in \code{y} = \code{TRUE} in \code{\link{lm}} function. To read more about this see the help documentation for
#' \code{\link{lm}} function or see the example below.
#'
#' @param ngrid the number of equally spaced points to discretize the (0,1) interval for computing the covariance function.
#'
#' @param gridpit logical. If \code{TRUE} (the default value), the parameter ngrid is ignored and (0,1) interval is divided
#' based on probability inverse transformed values obtained from the sample. If \code{FALSE}, the interval is divided into ngrid
#' equally spaced points for computing the covariance function.
#'
#' @param hessian logical. If \code{TRUE} the Fisher information matrix is estimated by the observed Hessian Matrix based on
#' the sample. If \code{FALSE} (the default value) the Fisher information matrix is estimated by the variance of the
#' observed score matrix.
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
testLMNormal = function(x, y, fit = NULL, ngrid = length(y), gridpit = FALSE, hessian = FALSE, method = 'cvm'){


  if( is.null(fit) ){

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

    # Apply linear model with normal assumption
    temp   <- applyLMNormal(x = x, y = y)

    # Extract score function, pit values, and MLE
    Score  <- temp$Score
    pit    <- temp$pit
    par    <- temp$par

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

    # Compute Cramer-von-Mises statistic
    cvm      <- getCvMStatistic(pit)

    # Compute pvalue
    pvalue  <- getpvalue(u = cvm, eigen = ev)
    res     <- list(Statistic = cvm, pvalue = pvalue)

    return(res)

  }else{

    if (!inherits(fit, 'lm')){
      stop('The fit must be \'lm\' object.')
    }

    if( !is.matrix(fit$x) ){
      stop('fit object must have the design matrix corresponding to the model. Consider setting x = TRUE in lm function to return x matrix.')
    }

    if( !is.vector(fit$y) ){
      stop('fit object must contain the response variable. Consider setting y = TRUE in lm function to return reponse variable.')
    }

    x <- fit$x
    y <- fit$y

    # Apply linear model with normal assumption
    temp   <- applyLMNormal(x = x, y = y)

    # Extract score function, pit values, and MLE
    Score  <- temp$Score
    pit    <- temp$pit
    par    <- temp$par

    # Compute Fisher information matrix
    if( hessian ){
      fisher <- lmFisherByHessian(x = x, y = y, theta = par)
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

      # Compute C-v-M statistic
      cvm      <- getCvMStatistic(pit)

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

      # Compute C-v-M statistic
      AD      <- getCvMStatistic(pit)

      # Compute pvalue
      pvalue  <- getpvalue(u = AD, eigen = ev)

      res     <- list(Statistic = AD, pvalue = pvalue)

      return(res)


    } else {

      # Calculate both cvm and ad statistics

      # 1. Do cvm calculation
      cvm        <- getCvMStatistic(pit)
      names(cvm) <- 'Cramer-von-Mises Statistic'

      # Compute Eigen values
      if( gridpit ){
        ev    <- getEigenValues(S = Score, FI = fisher, pit = pit, me = 'cvm')
      }else{
        ev    <- getEigenValues_manualGrid(S = Score, FI = fisher, pit = pit, M = ngrid, me = 'cvm')
      }

      # Compute pvalue
      cvm.pvalue  <- getpvalue(u = cvm, eigen = ev)
      names(cvm.pvalue) <- 'pvalue for Cramer-von-Mises test'


      # 2. Do AD calculations
      AD      <- getADStatistic(pit)
      names(AD) <- 'Anderson-Darling Statistic'

      # Compute Eigen values
      if( gridpit ){
        ev    <- getEigenValues(S = Score, FI = fisher, pit = pit, me = 'ad')
      }else{
        ev    <- getEigenValues_manualGrid(S = Score, FI = fisher, pit = pit, M = ngrid, me = 'ad')
      }

      # Calculate pvalue
      ad.pvalue  <- getpvalue(u = AD, eigen = ev)
      names(ad.pvalue) <- 'Anderson-Darling test'


      # Prepare a list to return both statistics and their approximate pvalue
      res     <- list(Statistics = c(cvm, AD), pvalue = c(cvm.pvalue, ad.pvalue) )
      return(res)

    }

  }

}
