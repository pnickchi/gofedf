#' Apply Goodness of Fit Test to the Residuals of a Generalized Linear Model with Gamma Link Function
#'
#'\code{testGLMGamma} is used to check the validity of Gamma assumption for the response variable when
#'fitting generalized linear model. Common link functions in \code{\link{glm}} can be used here.
#'
#' @param x is either a numeric vector or a design matrix. In the design matrix, rows indicate observations and columns
#' presents covariats.
#'
#' @param y is a vector of numeric values with the same number of observations or number of rows as x.
#'
#' @param l a character vector indicating the link function that should be used for Gamma family. Some common
#' link functions for Gamma family are 'log' and 'inverse'. For more details see \code{\link{make.link}} from stats
#' package in R.
#'
#' @param fit is an object of class \code{glm} and its default value is NULL. If a fit of class \code{glm} is provided,
#' the arguments \code{x}, \code{y}, and \code{l} will be ignored. We recommend using \code{\link{glm2}} function from
#' \code{\link{glm2}} package since it provides better convergence while optimizing the likelihood to estimate
#' coefficients of the model by IWLS method. It is required to return design matrix by \code{x} = \code{TRUE} in
#' \code{\link{glm}} or \code{\link{glm2}} function. For more information on how to do this, refer to the help
#' documentation for the \code{\link{glm}} or \code{\link{glm2}} function.
#'
#' @param start.value a numeric value or vector. This is the same as \code{start} argument in \code{\link{glm}} or
#' \code{\link{glm2}}. The value is a starting point in iteratively reweighted least squares (IRLS) algorithm for
#' estimating the MLE of coefficients in the model.
#'
#' @param control a list of parameters to control the fitting process in \code{glm} or \code{glm2} function.
#' For more details, see \code{\link{glm.control}}.
#'
#' @param method a character string indicating which goodness-of-fit statistic is to be computed. The default value is
#' 'cvm' for the Cramer-von-Mises statistic. Other options include 'ad' for the Anderson-Darling statistic, and 'both'
#' to compute both cvm and ad.
#'
#' @return A list of three containing the following components:
#' - Statistic: the value of goodness-of-fit statistic.
#' - pvalue: the approximate pvalue for the goodness-of-fit test based on empirical distribution function.
#' if method = 'cvm' or method = 'ad', it returns a numeric value for the statistic and pvalue. If method = 'both', it
#' returns a numeric vector with two elements and one for each statistic.
#' - converged: logical to indicate if the IWLS algorithm have converged or not.
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' n <- 50
#' p <- 5
#' x <- matrix( rnorm(n*p, mean = 10, sd = 0.1), nrow = n, ncol = p)
#' b <- runif(p)
#' e <- rgamma(n, shape = 3)
#' y <- exp(x %*% b) * e
#' testGLMGamma(x, y, l = 'log')
#' myfit <- glm(y ~ x, family = Gamma('log'), x = TRUE)
#' testGLMGamma(fit = myfit)
testGLMGamma = function(x, y, l, fit = NULL, start.value = NULL, control = NULL, method = 'cvm'){


  if( is.null(control) ){
    control <- glm.control(epsilon = 1e-8, maxit = 100, trace = F)
  }else{
    control <- control
  }

  if( is.null(fit) ){

    family  <- Gamma(link = l)
    n       <- length(y)
    temp    <- applyGLMGamma(x, y, fml = family, sv = start.value, ctl = control, fit.included = NULL)
    Score   <- temp$Score
    pit     <- temp$pit
    par     <- temp$par

    fisher  <- (n-1)*var(Score)/n

    if( method == 'cvm'){

      ev      <- getEigenValues(S = Score, FI = fisher, pit, me = 'cvm')
      cvm     <- getCvMStatistic(pit)
      pvalue  <- getpvalue(u = cvm, eigen = ev)

      res     <- list(Statistic = cvm, pvalue = pvalue, converged = temp$converged)
      return(res)

    } else if ( method == 'ad') {

      ev      <- getEigenValues(S = Score, FI = fisher, pit, me = 'ad')
      ad      <- getADStatistic(pit)
      pvalue  <- getpvalue(u = ad, eigen = ev)

      res     <- list(Statistic = U2, pvalue = pvalue, converged = temp$converged)
      return(res)

    } else {

      # Calculate both cvm and ad statisitcs

      # 1. Do cvm calculation
      cvm        <- getCvMStatistic(pit)
      names(cvm) <- 'Cramer-von-Mises Statistic'

      # Get Eigen values
      ev    <- getEigenValues(S = Score, FI = fisher, pit = pit, me = 'cvm')

      # Calculate pvalue
      cvm.pvalue  <- getpvalue(u = cvm, eigen = ev)
      names(cvm.pvalue) <- 'pvalue for Cramer-von-Mises test'


      # 2. Do ad calculations
      ad      <- getADStatistic(pit)
      names(ad) <- 'Anderson-Darling Statistic'

      # Get Eigen values
      ev    <- getEigenValues(S = Score, FI = fisher, pit = pit, me = 'ad')

      # Calculate pvalue
      ad.pvalue  <- getpvalue(u = ad, eigen = ev)
      names(ad.pvalue) <- 'Anderson-Darling test'


      # Prepare a list to return both statistics and their approximate pvalue
      res     <- list(Statistics = c(cvm, ad), pvalue = c(cvm.pvalue, ad.pvalue) )
      return(res)
    }


  }else{

    if (!inherits(fit, 'glm')){
      stop('The fit must be \'glm\' object returned by either glm or glm2 function.')
    }

    if( fit$family$family != 'Gamma' ){
      stop('The family in fit must be Gamma.')
    }

    if( !is.matrix(fit$x) ){
      stop('fit object must have the design matrix corresponding to the model. Consider setting x = TRUE in glm function to return x matrix.')
    }

    x       <- fit$x
    y       <- fit$y
    family  <- fit$family
    n       <- length(y)
    temp    <- applyGLMGamma(x, y, fml = family, sv = start.value, ctl = control, fit.included = fit)
    Score   <- temp$Score
    pit     <- temp$pit
    par     <- temp$par

    fisher  <- (n-1)*var(Score)/n


    if( method == 'cvm'){

      ev    <- getEigenValues(S = Score, FI = fisher, pit, me = 'cvm')
      U2      <- getCvMStatistic(pit)
      pvalue  <- getpvalue(u = U2, eigen = ev)

      res     <- list(Statistic = U2, pvalue = pvalue, converged = temp$converged)
      return(res)

    } else if ( method == 'ad') {

      ev      <- getEigenValues(S = Score, FI = fisher, pit, me = 'ad')
      U2      <- getCvMStatistic(pit)
      pvalue  <- getpvalue(u = U2, eigen = ev)

      res     <- list(Statistic = U2, pvalue = pvalue, converged = temp$converged)
      return(res)

    } else {

      # Calculate both cvm and ad statisitcs

      # 1. Do cvm calculation
      cvm        <- getCvMStatistic(pit)
      names(cvm) <- 'Cramer-von-Mises Statistic'

      # Get Eigen values
      ev    <- getEigenValues(S = Score, FI = fisher, pit = pit, me = 'cvm')

      # Calculate pvalue
      cvm.pvalue  <- getpvalue(u = cvm, eigen = ev)
      names(cvm.pvalue) <- 'pvalue for Cramer-von-Mises test'


      # 2. Do ad calculations
      ad      <- getADStatistic(pit)
      names(ad) <- 'Anderson-Darling Statistic'

      # Get Eigen values
      ev    <- getEigenValues(S = Score, FI = fisher, pit = pit, me = 'ad')

      # Calculate pvalue
      ad.pvalue  <- getpvalue(u = ad, eigen = ev)
      names(ad.pvalue) <- 'Anderson-Darling test'


      # Prepare a list to return both statistics and their approximate pvalue
      res     <- list(Statistics = c(cvm, ad), pvalue = c(cvm.pvalue, ad.pvalue), converged = temp$converged )
      return(res)

    }

  }


}
