applyGLMGamma = function(x, y, fml, sv, ctl, fit.included){

   # Check if we need to fit a GLM or user already provided the fit
   if( is.null(fit.included) ){

     # Fit a generalized linear model by glm2 package
     fit <- glm2::glm2(formula = y ~ x,
                       family = fml,
                       x = TRUE,
                       control = ctl,
                       start   = sv,
                       na.action = na.omit)

   }else{
     fit <- fit.included
   }


   ###############################################################################################

   #
   # This section calculates the MLE of parameters in the model.
   # There are p coefficients in beta vector.
   # We need to estimate the MLE of shape parameter as well.
   #

   # coef function from stats package extracts MLE estimate of beta parameter in the model.
   mle.coef                 <- coef(fit)

   #
   mle.alpha                <- MASS::gamma.shape(fit, it.lim = 30, eps.max = 1e-9)[1]$alpha

   # Build a vector for MLE estiamtes, p parameters for beta and last parameter for mle of shape parameter
   par                      <- c(mle.coef, mle.alpha)
   names(par)[length(par)]  <- 'shape'

   ###############################################################################################



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
