% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testGLMGamma.R
\name{testGLMGamma}
\alias{testGLMGamma}
\title{Apply Goodness of Fit Test to the Residuals of a Generalized Linear Model with Gamma Link Function}
\usage{
testGLMGamma(
  x,
  y,
  fit = NULL,
  l = "log",
  hessian = FALSE,
  start.value = NULL,
  control = NULL,
  method = "cvm"
)
}
\arguments{
\item{x}{is either a numeric vector or a design matrix. In the design matrix, rows indicate observations and columns
presents covariats.}

\item{y}{is a vector of numeric values with the same number of observations or number of rows as x.}

\item{fit}{is an object of class \code{glm} and its default value is NULL. If a fit of class \code{glm} is provided,
the arguments \code{x}, \code{y}, and \code{l} will be ignored. We recommend using \code{\link{glm2}} function from
\code{\link{glm2}} package since it provides better convergence while optimizing the likelihood to estimate
coefficients of the model by IWLS method. It is required to return design matrix by \code{x} = \code{TRUE} in
\code{\link{glm}} or \code{\link{glm2}} function. For more information on how to do this, refer to the help
documentation for the \code{\link{glm}} or \code{\link{glm2}} function.}

\item{l}{a character vector indicating the link function that should be used for Gamma family. Some common
link functions for Gamma family are 'log' and 'inverse'. For more details see \code{\link{make.link}} from stats
package in R.}

\item{hessian}{logical. If \code{TRUE} the Fisher information matrix is estimated by the observed Hessian Matrix based on
the sample. If \code{FALSE} (the default value) the Fisher information matrix is estimated by the variance of the
observed score matrix.}

\item{start.value}{a numeric value or vector. This is the same as \code{start} argument in \code{\link{glm}} or
\code{\link{glm2}}. The value is a starting point in iteratively reweighted least squares (IRLS) algorithm for
estimating the MLE of coefficients in the model.}

\item{control}{a list of parameters to control the fitting process in \code{glm} or \code{glm2} function.
For more details, see \code{\link{glm.control}}.}

\item{method}{a character string indicating which goodness-of-fit statistic is to be computed. The default value is
'cvm' for the Cramer-von-Mises statistic. Other options include 'ad' for the Anderson-Darling statistic, and 'both'
to compute both cvm and ad.}
}
\value{
A list of three containing the following components:
\itemize{
\item Statistic: the value of goodness-of-fit statistic.
\item p-value: the approximate p-value for the goodness-of-fit test based on empirical distribution function.
if method = 'cvm' or method = 'ad', it returns a numeric value for the statistic and p-value. If method = 'both', it
returns a numeric vector with two elements and one for each statistic.
\item converged: logical to indicate if the IWLS algorithm have converged or not.
}
}
\description{
\code{testGLMGamma} is used to check the validity of Gamma assumption for the response variable when
fitting generalized linear model. Common link functions in \code{\link{glm}} can be used here.
}
\examples{
set.seed(123)
n <- 50
p <- 5
x <- matrix( rnorm(n*p, mean = 10, sd = 0.1), nrow = n, ncol = p)
b <- runif(p)
e <- rgamma(n, shape = 3)
y <- exp(x \%*\% b) * e
testGLMGamma(x, y, l = 'log')
myfit <- glm(y ~ x, family = Gamma('log'), x = TRUE, y = TRUE)
testGLMGamma(fit = myfit)
}
