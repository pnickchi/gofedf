#' Apply GLM to sample with the assumption of Gamma response and compute required components for the test.
#'
#' @description Compute maximum likelihood estimates of the coefficients and shape parameter in Gamma distribution, Score function
#' evaluated at the sample, and probability inverse transformed (PIT) values of sample.
#'
#' @param x a matrix of explanatory variables.
#'
#' @param y a numeric vector of response values.
#'
#' @param fml a character vector to define family function in GLM.
#'
#' @param sv starting value for computing MLE of coefficients.
#'
#' @param ctl a list of parameters to control the fitting process in \code{glm} or \code{glm2} function.
#'
#' @param fit.included logical to indicate if the fit is included as an object of class \code{glm}.
#'
#' @return a list with three elements.
#'
#' @noRd
applyGLMGamma = function(x, y, fml, sv, ctl, fit.included){

  #
  # If the fit is not included: fit a model with glm2
  #
  if( is.null(fit.included) ){

    # Check if the starting values are included
    if( is.null(sv) ){
      fit <- glm2::glm2(formula = y ~ x, family = fml, x = TRUE, control = ctl, na.action = na.omit)
    }

    if( !is.null(sv) ){
      fit <- glm2::glm2(formula = y ~ x, family = fml, x = TRUE, control = ctl, start = sv, na.action = na.omit)
    }

   }


   if( !is.null(fit.included)){
    fit <- fit.included
   }


   #
   # This section calculates the MLE of parameters in the model.
   # There are p coefficients in beta vector.
   # We need to estimate the MLE of shape parameter as well.
   # In total, there is (p+1) parameters in the model.

   # coef function from stats package extracts MLE estimate of beta parameter in the model.
   mle.coef                 <- coef(fit)

   # Use gamma.shape function from MASS package to compute mle of shape parameter in Gamma distribution
   mle.alpha                <- MASS::gamma.shape(fit, it.lim = 30, eps.max = 1e-9)[1]$alpha

   # Build a vector for MLE estimates, p parameters for beta and last parameter for mle of shape parameter
   par                      <- c(mle.coef, mle.alpha)
   names(par)[length(par)]  <- 'shape'


   # Extract family function and design matrix (X) from fit object and calculates the linear predictor.
   fm              <- family(fit)
   X               <- fit$x
   linearPredictor <- X %*% mle.coef

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
