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
#' @param fit is an object of class \code{glm} and its default value is NULL. If a fit of class \code{glm} is provided,
#' the arguments \code{x}, \code{y}, and \code{l} will be ignored. We recommend using \code{\link{glm2}} function from
#' \code{\link{glm2}} package since it provides better convergence while optimizing the likelihood to estimate
#' coefficients of the model by IWLS method. It is required to return design matrix by \code{x} = \code{TRUE} in
#' \code{\link{glm}} or \code{\link{glm2}} function. For more information on how to do this, refer to the help
#' documentation for the \code{\link{glm}} or \code{\link{glm2}} function.
#'
#' @param l a character vector indicating the link function that should be used for Gamma family. Some common
#' link functions for Gamma family are 'log' and 'inverse'. For more details see \code{\link{make.link}} from stats
#' package in R.
#'
#' @param hessian logical. If \code{TRUE} the Fisher information matrix is estimated by the observed Hessian Matrix based on
#' the sample. If \code{FALSE} (the default value) the Fisher information matrix is estimated by the variance of the
#' observed score matrix.
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
#' - p-value: the approximate p-value for the goodness-of-fit test based on empirical distribution function.
#' if method = 'cvm' or method = 'ad', it returns a numeric value for the statistic and p-value. If method = 'both', it
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
#' myfit <- glm(y ~ x, family = Gamma('log'), x = TRUE, y = TRUE)
#' testGLMGamma(fit = myfit)
testGLMGamma = function(x, y, fit = NULL, l = 'log', hessian = FALSE, start.value = NULL, control = NULL, method = 'cvm'){

   if( is.null(fit) ){

     # Check if the link is valid
     if( !(l %in% c('inverse','identity','log')) ){
       stop('The link for Gamma must be inverse, identity, or log.')
     }

     # Assign control vector for glm2 function
     if( is.null(control) ){
       ctl <- glm.control(epsilon = 1e-8, maxit = 100, trace = F)
     }else{
       ctl <- control
     }

     # Check if the starting value is provided
     if( is.null(start.value) ){

       # Fit a GLM-Gamma to the data to compute maximum likelihood estimation of coefficients
       fitobj <- glm2::glm2(formula = y ~ x, family = Gamma(link = l), x = TRUE, y = TRUE, control = ctl, na.action = na.omit)

     }else{

       if( length(start.value) != (ncol(x)+1) ){
         stop('The lenght of starting.value does not match the number of columns in x.')
       }

       # Fit a GLM-Gamma to the data to compute maximum likelihood estimation of coefficients
       fitobj <- glm2::glm2(formula = y ~ x, family = Gamma(link = l), x = TRUE, y = TRUE, control = ctl, na.action = na.omit, start = start.value)
     }

     # Get the sample size
     n       <- length(fitobj$y)

   }

   if( !is.null(fit) ){

     if (!inherits(fit, 'glm')){
       stop('The fit must be \'glm\' object returned by either glm or glm2 function.')
     }

     if( fit$family$family != 'Gamma' ){
       stop('The family must be Gamma.')
     }

     if( !(fit$family$link %in% c('inverse','identity','log')) ){
       stop('The link for Gamma must be inverse, identity, or log.')
     }

     if( !is.matrix(fit$x) ){
       stop('fit must contain the design matrix. Consider setting x = TRUE in glm or glm2 function to return x matrix.')
     }

     if( !is.vector(fit$y) ){
       stop('fit must contain the response variable. Consider setting y = TRUE in glm or glm2 function to return response variable.')
     }

     fitobj <- fit
   }

   # Apply GLM Gamma to compute score, MLE of parameters, and pit values
   temp    <- applyGLMGamma(fit = fitobj)
   Score   <- temp$Score
   pit     <- temp$pit
   par     <- temp$par

   # Boolean to check if the IWLS algorithm have converged
   converged <- temp$converged

   if( !converged ){
     message('The IWLS iterative algorithm did not converge.')
   }

   # Compute Fisher information matrix
   if(hessian){
     fisher  <- glmgammaFisherByHessian(fit = fitobj, mle_shape = par['shape'])
   }else{
     fisher  <- (n-1)*var(Score)/n
   }

   if( method == 'cvm'){

     # Compute Eigen values
     ev      <- getEigenValues(S = Score, FI = fisher, pit, me = 'cvm')

     # Compute cvm statistic
     cvm     <- getCvMStatistic(pit)

     # Compute p-value
     pvalue  <- getpvalue(u = cvm, eigen = ev)

     res     <- list(Statistic = cvm, pvalue = pvalue, converged = converged)
     return(res)

    } else if ( method == 'ad') {

     # Compute Eigen values
     ev      <- getEigenValues(S = Score, FI = fisher, pit, me = 'ad')

     # Compute ad statistics
     ad      <- getADStatistic(pit)

     # Compute p-value
     pvalue  <- getpvalue(u = ad, eigen = ev)

     res     <- list(Statistic = ad, pvalue = pvalue, converged = converged)
     return(res)

    } else {

      # Calculate both cvm and ad statisitcs

      # Compute cvm statistic
      cvm        <- getCvMStatistic(pit)
      names(cvm) <- 'Cramer-von-Mises Statistic'

      # Compute Eigen values
      ev    <- getEigenValues(S = Score, FI = fisher, pit = pit, me = 'cvm')

      # Compute p-value
      cvm.pvalue  <- getpvalue(u = cvm, eigen = ev)
      names(cvm.pvalue) <- 'pvalue for Cramer-von-Mises test'


      # 2. Compute ad calculations
      ad      <- getADStatistic(pit)
      names(ad) <- 'Anderson-Darling Statistic'

      # Compute Eigen values
      ev    <- getEigenValues(S = Score, FI = fisher, pit = pit, me = 'ad')

      # Compute p-value
      ad.pvalue  <- getpvalue(u = ad, eigen = ev)
      names(ad.pvalue) <- 'Anderson-Darling test'


      # Prepare a list to return both statistics and their approximate pvalue
      res     <- list(Statistics = c(cvm, ad), pvalue = c(cvm.pvalue, ad.pvalue) )
      return(res)
    }

}
