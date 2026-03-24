#'Apply Goodness of Fit Test to the Residuals of a Generalized Linear Model with
#'Gamma Link Function
#'
#'\code{testGLMGamma} is used to check the validity of Gamma assumption for the
#'response variable when fitting generalized linear model. Common link functions
#'in \code{\link{glm}} can be used here.
#'
#'@param x is either a numeric vector or a design matrix. In the design matrix,
#'  rows indicate observations and columns presents covariats.
#'
#'@param y is a vector of numeric values with the same number of observations or
#'  number of rows as x.
#'
#'@param fit is an object of class \code{glm} and its default value is NULL. If
#'  a fit of class \code{glm} is provided, the arguments \code{x}, \code{y}, and
#'  \code{l} will be ignored. We recommend using \code{\link[glm2]{glm2}} function
#'  from \code{\link[glm2]{glm2}} package since it provides better convergence while
#'  optimizing the likelihood to estimate coefficients of the model by IWLS
#'  method. It is required to return design matrix by \code{x} = \code{TRUE} in
#'  \code{\link{glm}} or \code{\link[glm2]{glm2}} function. For more information on
#'  how to do this, refer to the help documentation for the \code{\link{glm}} or
#'  \code{\link[glm2]{glm2}} function.
#'
#'@param l a character vector indicating the link function that should be used
#'  for Gamma family. Acceptable link functions for Gamma family are inverse,
#'  identity and log. For more details see \code{\link{Gamma}} from stats
#'  package.
#'
#'@param discretize If \code{TRUE}, the covariance function of \eqn{W_{n}(u)}
#'  process is evaluated at some data points (see \code{ngrid} and
#'  \code{gridpit}), and the integral equation is replaced by a matrix equation.
#'  If \code{FALSE} (the default value), the covariance function is first
#'  estimated, and then the integral equation is solved to find the eigenvalues.
#'  The results of our simulations recommend using the estimated covariance for
#'  solving the integral equation. The parameters \code{ngrid}, \code{gridpit},
#'  and \code{hessian} are only relevant when \code{discretize = TRUE}.
#'
#'@param ngrid the number of equally spaced points to discretize the (0,1)
#'  interval for computing the covariance function.
#'
#' @param gridpit logical. If \code{TRUE} (the default value), the parameter
#'   ngrid is ignored and (0,1) interval is divided based on probability
#'   integral transforms or PITs obtained from the sample. If \code{FALSE}, the
#'   interval is divided into ngrid equally spaced points for computing the
#'   covariance function.
#'
#'@param hessian logical. If \code{TRUE} the Fisher information matrix is
#'  estimated by the observed Hessian Matrix based on the sample. If
#'  \code{FALSE} (the default value) the Fisher information matrix is estimated
#'  by the variance of the observed score matrix.
#'
#'@param start.value a numeric value or vector. This is the same as \code{start}
#'  argument in \code{\link{glm}} or \code{\link[glm2]{glm2}}. The value is a starting
#'  point in iteratively reweighted least squares (IRLS) algorithm for
#'  estimating the MLE of coefficients in the model.
#'
#'@param control a list of parameters to control the fitting process in
#'  \code{glm} or \code{glm2} function. For more details, see
#'  \code{\link{glm.control}}.
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
#'@export
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
#' # Example for custom weight function
#' w_cvm <- function(u) rep(1, length(u))
#' testGLMGamma(fit = myfit, method = 'user', weight_function = w_cvm)
#'
testGLMGamma = function(x, y, fit = NULL, l = 'log', discretize = FALSE, ngrid = length(y), gridpit = TRUE, hessian = FALSE, start.value = NULL, control = NULL, method = 'cvm', weight_function = NULL){

   if( is.null(fit) ){

     # Make sure all observations in response are positive
     if( any(y <= 0) ){
       stop('y values must be positive for Gamma distribution.')
     }

     # Check if the link is valid
     if( !(l %in% c('inverse','identity','log')) ){
       stop('The link for Gamma must be either inverse, identity, or log.')
     }

     # Make sure discretize is logical type
     if( !is.logical(discretize) ){
       stop('discretize must be either TRUE or FALSE.')
     }

     # Make sure ngrid is integer type
     if( !(ngrid %% 1 == 0) ){
       stop('ngrid must be an integer number.')
     }

     # Make sure discretize and hessian are logical type
     if( !is.logical(gridpit) ){
       stop('gridpit must be either TRUE or FALSE.')
     }

     if( !is.logical(hessian) ){
       stop('hessian must be either TRUE or FALSE.')
     }

     if( !(method %in% c('cvm','ad','both','user')) ){
       stop('method must be either cvm, ad, both, or user.')
     }

     if( method == 'user' & is.null(weight_function) ){
       stop('method is set to user but weight_function is NULL. You have to pass a weight function if you use user method.')
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

     mu <- fitobj$family$linkfun(fitobj$linear.predictors)

     if( any(mu < 0) ){
       stop('The fitted mean response is negative.')
     }

   }

   if( !is.null(fit) ){

     if (!inherits(fit, 'glm')){
       stop('The fit must be \'glm\' object returned by either glm function (from stats package) or glm2 function (from glm2 package).')
     }

     if( !(method %in% c('cvm','ad','both','user')) ){
       stop('method must be either cvm, ad, both, user.')
     }

     if( fit$family$family != 'Gamma' ){
       stop('The family must be Gamma.')
     }

     if( !(fit$family$link %in% c('inverse','identity','log')) ){
       stop('The link for Gamma must be inverse, identity, or log.')
     }

     if( !is.matrix(fit$x) | !is.vector(fit$y) ){
       stop('fit must contain the design matrix and the response variable.  \n Consider setting x = TRUE and y = TRUE in glm or glm2 function to return both.')
     }

     fitobj <- fit

     mu <- fitobj$family$linkfun(fitobj$linear.predictors)

     if( any(mu < 0) ){
       stop('The fitted mean response is negative.')
     }

   }

   # Apply GLM Gamma to compute score, MLE of parameters, and pit values
   temp    <- applyGLMGamma(fit = fitobj)
   Score   <- temp$Score
   pit     <- temp$pit
   par     <- temp$par

   # Get the sample size
   n <- nrow(Score)

   # Boolean to check if the IWLS algorithm have converged
   converged <- temp$converged

   if( !converged ){
     message('The IWLS iterative algorithm did not converge. \n Consider increasing maxit or decreasing epsilon in glm.control() function.')
   }


   #
   # Use the estimated covariance function and solving the integral equation analytically
   #
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

     if( method == 'user' ){

       # Compute P matrix
       P <- computeMatrix(n, Score, method = method, w_function = weight_function)

       # Adjust for the number of estimated parameters
       P <- P / (n-ncol(Score)-1)

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
   # Use the estimated covariance function and turning integral equation into a matrix equation
   #

   # Compute Fisher information matrix
   if(hessian){
     fisher  <- glmgammaFisherByHessian(fit = fitobj, mle_shape = par['shape'])
   }else{
     fisher  <- (n-1)*var(Score)/n
   }

   if( method == 'cvm'){

     # Compute Eigen values
     if( gridpit ){
       ev    <- getEigenValues(S = Score, FI = fisher, pit, me = method)
     }else{
       ev    <- getEigenValues_manualGrid(S = Score, FI = fisher, pit, M = ngrid, me = method)
     }

     # Compute cvm statistic
     cvm     <- getCvMStatistic(pit)

     # Compute p-value
     pvalue  <- getpvalue(u = cvm, eigen = ev)

     res     <- list(Statistic = cvm, pvalue = pvalue, converged = converged)
     return(res)

    } else if ( method == 'ad') {

     # Compute Eigen values
     if( gridpit ){
       ev    <- getEigenValues(S = Score, FI = fisher, pit, me = method)
     }else{
       ev    <- getEigenValues_manualGrid(S = Score, FI = fisher, pit, M = ngrid, me = method)
     }

     # Compute ad statistics
     ad      <- getADStatistic(pit)

     # Compute p-value
     pvalue  <- getpvalue(u = ad, eigen = ev)

     res     <- list(Statistic = ad, pvalue = pvalue, converged = converged)
     return(res)

    }else if ( method == 'both' ){

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

    }else{

      # Compute Eigen values
      if( gridpit ){
        ev    <- getEigenValues(S = Score, FI = fisher, pit, me = method, w_function = weight_function)
      }else{
        ev    <- getEigenValues_manualGrid(S = Score, FI = fisher, pit, M = ngrid, me = method, w_function = weight_function)
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
