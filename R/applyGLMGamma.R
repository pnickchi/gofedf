#' Apply GLM to sample with the assumption of Gamma response and compute required components for the test.
#'
#' @description Compute maximum likelihood estimates of the coefficients and shape parameter in Gamma distribution, Score function
#' evaluated at the sample, and probability inverse transformed (PIT) values of sample.
#'
#' @param fit is an object of class \code{glm} and its default value is NULL. If a fit of class \code{glm} is provided,
#' the arguments \code{x}, \code{y}, and \code{l} will be ignored. We recommend using \code{\link{glm2}} function from
#' \code{\link{glm2}} package since it provides better convergence while optimizing the likelihood to estimate
#' coefficients of the model by IWLS method. It is required to return design matrix by \code{x} = \code{TRUE} in
#' \code{\link{glm}} or \code{\link{glm2}} function. For more information on how to do this, refer to the help
#' documentation for the \code{\link{glm}} or \code{\link{glm2}} function.
#'
#' @return a list with three elements.
#'
#' @noRd
applyGLMGamma = function(fit){

   # coef function from stats package extracts MLE estimate of beta parameter in the model.
   mle.coef                 <- coef(fit)

   # Use gamma.shape function from MASS package to compute MLE of shape parameter in Gamma distribution
   # Note: in our testing and simluation, we did not see any improvement in the estimate by increasing
   # it.lim or reducing eps.max.
   mle.alpha                <- MASS::gamma.shape(fit, it.lim = 30, eps.max = 1e-9)[1]$alpha

   # Build a vector for MLE estimates, p parameters for beta and last parameter for mle of shape parameter
   par                      <- c(mle.coef, mle.alpha)
   names(par)[length(par)]  <- 'shape'


   # Extract family function, design matrix (X) from fit object, compute the linear predictor and extracy response variable.
   fm              <- family(fit)
   X               <- fit$x
   linearPredictor <- X %*% mle.coef
   y               <- fit$y

   # Check if we need to add offset to the linear predictor.
   if( !is.null(fit$offset) ){
     linearPredictor <- linearPredictor + fit$offset
   }

   # Calculate miohat, estimated value of mean function in GLM
   miohat          <- fm$linkinv(linearPredictor)

   # Calculate the score function
   # Use family function to get link function and derivative of the inverse-link function
   S1    <- 1 - digamma(mle.alpha) + log(mle.alpha) + log(y/miohat) - (y/miohat)
   S2    <- ( fm$mu.eta(linearPredictor) / miohat ) * mle.alpha * ( (y/miohat) - 1 )
   S2    <- as.vector(S2) * X
   S     <- cbind(S1,S2)

   # Calculate probability inverse transfer of data
   pit <- pgamma(q = as.numeric(y), shape = mle.alpha, scale = miohat/mle.alpha)

   # Define the list to return
   res <- list(Score = S, pit = pit, par = par, converged = fit$converged)

   # Return results
   return(res)

}
