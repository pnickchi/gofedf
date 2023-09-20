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

   # Compute MLE of parameters (coefficients in the model + shape parameter in Gamma)
   par <- glmMLE(fit)

   # Compute the Score and PIT for the case of GLM with Gamma
   temp <- glmScorePIT(fit, theta = par['shape'])
   S    <- temp$Score
   pit  <- temp$pit

   # Define the list to return
   res <- list(Score = S, pit = pit, par = par, converged = fit$converged)

   # Return results
   return(res)

}
