README
================
2023-07-25

# gofedf

<!-- badges: start -->

[![CRANstatus](https://www.r-pkg.org/badges/version/gofedf)](https://cran.r-project.org/package=gofedf)
[![GitHubRelease](https://img.shields.io/github/release/pnickchi/gofedf?style=flat)](https://github.com/pnickchi/gofedf/releases)
[![Github All
Releases](https://img.shields.io/github/downloads/pnickchi/gofedf/total.svg?style=flat)](https://github.com/pnickchi/gofedf)
<!-- badges: end -->

Authors:

- [Richard Lockhart](http://www.sfu.ca/~lockhart/), <lockhart@sfu.ca>
- [Payman Nickchi](https://github.com/pnickchi), <pnickchi@sfu.ca>
  (Maintainer)

The `gofedf` package provides computational tools to apply
goodness-of-fit tests based on empirical distribution function theory.
The package offers functions and routines to test the hypothesis that a
univariate sample follows a distribution based on the empirical
distribution function. The theory is founded on reducing the problem to
a stochastic process and computing its covariance function. An
approximate p-value is computed using the Farebrother method based on
the limiting distribution of the statistic. Users can run the test by
calculating either the Cramer-von Mises or Anderson-Darling statistic.
The covariance function of the stochastic process relies on specific
characteristics of the assumed model. Notably, knowledge of the Fisher
information matrix and the partial derivatives of the cumulative
distribution function is crucial for computing the covariance function.
However, obtaining these quantities can be computationally intensive or
challenging in general likelihood models. To overcome this limitation,
we propose an alternative method for estimating the covariance function
of the stochastic process directly from the sample data. The package
provides tools for this estimation and for testing if a sample comes
from any general likelihood model. In summary, the package can be used
to apply a goodness-of-fit test in any of the following settings:

- Validate if the assumptions about the response variable in a
  generalized linear model (with any link function) are satisfied. The
  current version only checks for the Gamma distribution.

- Validate if the normality assumptions in a linear model are satisfied.

- Apply a goodness-of-fit test to examine if a set of bivariate samples
  follows a Normal, Gamma, or Exponential distribution.

- Formal model evaluation in a general likelihood model. In this case,
  probability inverse transformed (PIT) values of the sample and the
  score function (if any parameter estimation is involved) are required.
  See the example for more details.

## Installation

You can install the released version of `gofedf` from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages('gofedf')
```

or the development version from [GitHub](https://github.com/pnickchi)
page with:

``` r
# install.packages("devtools")
devtools::install_github('pnickchi/gofedf')
```

## Examples

In this section, we will review the package’s primary functions and
provide illustrative examples to demonstrate their usage. The first two
examples relate to applying the goodness-of-fit (GOF) test for i.i.d.
samples, while the last two focus on linear models and generalized
linear models. The final example showcases the most important feature of
the package, allowing you to apply goodness-of-fit tests based on
empirical distribution function EDF for a general likelihood model.

### 1. Bivariate Normal distribution

The first example illustrates the GOF test for an i.i.d. sample from a
Normal distribution. The main function is `testNormal`. At the minimum,
it requires a numeric vector as input. By default, it uses probability
inverse-transformed values to compute the stochastic process and its
covariance function later. You can change this behavior by setting
`gridpit = FALSE` and assigning a positive value for `ngrid`. The
default value for `ngrid` is the same as the number of observations,
`n`. This means the (0,1) interval is divided into `n` equally spaced
data points to compute the stochastic process and the covariance
function. Additionally, the Fisher information matrix, by default, is
estimated by the variance of the score function. To change this, you can
set `hessian=TRUE` to estimate the Fisher information matrix using the
Hessian instead. Finally, method is a string that defines the statistic
to compute. Possible values are `cvm` for Cramer-von-Mises, `ad` for
Anderson-Darling, and both to compute both.

``` r
# Reproducible example
set.seed(123)

# Randomly generate some data from Normal distribution
n <- 50
x <- rnorm(n)

# Test if the data follows a Normal distribution by calculating the Cramer-von Mises statistic and approximate p-value of the test.
testNormal(x = x, method = 'cvm')
```

    ## $Statistic
    ## [1] 0.03781322
    ## 
    ## $pvalue
    ## [1] 0.6646717

``` r
# Test if the data follows a Normal distribution by calculating the Anderson-Darling statistic and approximate p-value of the test.
testNormal(x = x, method = 'ad')
```

    ## $Statistic
    ## [1] 0.2179704
    ## 
    ## $pvalue
    ## [1] 0.7833757

``` r
# Generate some random sample from a non Normal distribution.
x <- rgamma(n, shape = 3)
testNormal(x = x, method = 'cvm')
```

    ## $Statistic
    ## [1] 0.2141872
    ## 
    ## $pvalue
    ## [1] 0.002816061

### 2. Bivariate Gamma distribution

The second example illustrates the GOF test for an i.i.d. sample from a
Gamma distribution. The main function is `testGamma` and the arguments
remain the same as Normal case.

``` r
# Reproducible example
set.seed(123)

# Randomly generate some data
n <- 50
x <- rgamma(n, shape = 3)

# Test if the data follows a Gamma distribution, calculate Cramer-von Mises statistic and approximate p-value
testGamma(x = x, method = 'cvm')
```

    ## $Statistic
    ## [1] 0.0549759
    ## 
    ## $pvalue
    ## [1] 0.3318735

``` r
# Generate some random sample from a distribution that is not Gamma
x <- runif(n)
testNormal(x = x, method = 'cvm')
```

    ## $Statistic
    ## [1] 0.07085577
    ## 
    ## $pvalue
    ## [1] 0.1464319

### 3. Linear model with Normal error terms

In this example, we illustrate how to apply GOF test to verify the
assumptions of a linear model. The main function is `testLMNormal`. At
the minimum, a numeric vector of response variable, `y`, and a vector/
matrix of explanatory variables, `x`, are required. Conveniently, the
function can take an object of class “lm” and directly applies the
goodness-of-fit test. In this case, there is no need to pass `x` and
`y`. Note that if you decide to use this feature, you need to explicitly
ask `lm` function to return the design matrix and response variable by
passing `x=TRUE` and `y=TRUE` (as shown in the example below). The other
arguments of the function are the same as previous examples.

``` r
# Reproducible example
set.seed(123)

# Create a set of explanatory variables and a response variable according to a linear model

# Sample size
n <- 50

# Number of explanatory variables
p <- 5

# Generate some coefficients
b <- runif(p)

# Simulate random explanatory variables
X <- matrix( runif(n*p), nrow = n, ncol = p)

# Generate some error terms from Normal distribution
e <- rnorm(n)

# Generate response variable according to the linear model
y <- X %*% b + e

# Test if the residuals of the model follows a Normal distribution, calculate Cramer-von Mises statistic and approximate p-value
testLMNormal(x = X, y)
```

    ## $Statistic
    ## [1] 0.02437418
    ## 
    ## $pvalue
    ## [1] 0.8740983

``` r
# Or alternatively just pass 'lm.fit' object directly instead:
lm.fit <- lm(y ~ X - 1, x = TRUE, y = TRUE)
testLMNormal(fit = lm.fit)
```

    ## $Statistic
    ## [1] 0.02437418
    ## 
    ## $pvalue
    ## [1] 0.8740983

### 4. Gamma GLM with any link function

In this example, we illustrate how to apply the GOF test to verify if
the response variable in a generalized linear model with any link
function follows a Gamma distribution. The main function for this is
`testGLMGamma`. At the minimum, you need a numeric vector of the
response variable, `y`, a vector/matrix of explanatory variables, `x`,
and a link function for the Gamma family. Conveniently, the function can
take an object of class `glm` and directly apply the goodness-of-fit
test. You can use `glm` or `glm2` function from `glm2` pacakge. We
recommend using the `glm2` function from the `glm2` package as it
provides better estimates for the coefficients and avoids convergence
issues in the optimization process. In either of these cases, there is
no need to pass `x` and `y`. However, if you decide to use this feature,
you must explicitly ask the `glm` or `glm2` function to return the
design matrix and response variable by passing `x=TRUE` and `y=TRUE` (as
shown in the example below). Additionally, you can pass a starting
value, `start.value`, to be used as the initial value for MLE estimation
of the coefficients. The function also offers a list of parameters to
control the fitting process in `glm` or `glm2` functions. The other
arguments of the function remain consistent with previous examples.

``` r
# Reproducible example
set.seed(123)

# Create a set of explanatory variables and a response variable according to a generalized linear model.

# Sample size
n <- 50

# Number of explanatory variables
p <- 5

# Simulate random explanatory variables
X <- matrix( rnorm(n*p, mean = 10, sd = 0.1), nrow = n, ncol = p)

# Generate some coefficients
b <- runif(p)

# Generate some error terms from Gamma distribution
e <- rgamma(n, shape = 3)

# Generate response variable according to the generalized linear model (log link function)
y <- exp(X %*% b) * e

# Test if the Gamma assumptions of the response variable holds by calculating the Cramer-von Mises statistic and approximate p-value
testGLMGamma(x=X, y, l = 'log', method = 'cvm')
```

    ## $Statistic
    ## [1] 0.0870493
    ## 
    ## $pvalue
    ## [1] 0.1139862
    ## 
    ## $converged
    ## [1] TRUE

``` r
# Or alternatively just pass 'glm.fit' object directly instead:
glm.fit <- glm2::glm2(y ~ X, family = Gamma(link = 'log'), x = TRUE, y = TRUE)
testGLMGamma(fit = glm.fit, l = 'log')
```

    ## $Statistic
    ## [1] 0.0870493
    ## 
    ## $pvalue
    ## [1] 0.1139862
    ## 
    ## $converged
    ## [1] TRUE

### 5. General likelihood model

One of the most important features of the package is to provide
computational tools to apply the goodness-of-fit test based on empirical
distribution functions for any general likelihood model. We provided
tools to apply the test for Normal, Gamma, verify the assumptions in a
linear model and generalized linear model. But this additional feature
allows you to test if the sample come from any general likelihood model.
For example, consider you have a sample of size $n$, such as
$X_1, X_2, \ldots, X_n$, from a model with CDF of $F(X;\theta)$ where
$\theta$ contains $p$ parameters. Before running the test, at the
minimum you need the followings:

1)  A numeric vector of observations, `x`.

2)  The probability inverse transformed or PIT values of the sample
    which ought to be a numeric vector with the same size as `x` and
    with elements $F^{-1}(X_{i};\theta)$.

If $\theta$ is unknown, you also need to provide score function. This
needs to be a matrix with $n$ rows and $p$ columns where each row
measures the score of each observation. Note that the values are
computed as
$S(X_{i};\theta) = \frac{\partial}{\partial \theta} \log(f(X_{i};\theta))$
where $f(X_{i};\theta)$ is the probability density function. For sure,
if $\theta$ is not known, this means you need to compute the MLE of
$\theta$ to obtain item (ii) and if needed item (iii). The main function
to apply the GOF test in this case is `testYourModel`. The `precision`
argument sets the precision needed to check if the col sums of score
matrix are close enough to zero. The other arguments of the function
remain consistent with previous examples.

In the following example, we demonstrate how to apply the
goodness-of-fit test to check if a sample follows an Inverse Gaussian
distribution, where the shape parameter depends on some weights. First,
we generate data from an Inverse Gaussian distribution. For illustrative
purposes, we include functions to compute the Maximum Likelihood
Estimation (MLE) and score function for the sample. In the following
chunck of code, `IG_scorefunc` is a function that returns the score for
each observation, and `IG_pitfunc` is a function that provides a vector
of Probability Inverse Transformed (PIT) values. Additionally,
`IG_mlefunc` calculates the MLE of the mean and shape parameter. Second
we calculate score, PIT and MLE of parameters. Finally we call
`testYourModel` function to apply the test.

``` r
# Example: Inverse Gaussian (IG) distribution with weights

# Reproducible example
set.seed(123)

# Set the sample size
n <- 50

# Assign weights
covariates <- rep(1.5,n)

# Set mean and shape parameters for IG distribution.
mio        <- 2
lambda     <- 2

# Generate a random sample from IG distribution with weighted shape.
y <- statmod::rinvgauss(n, mean = mio, shape = lambda * covariates)

# Compute MLE of parameters, score matrix, and pit values.
theta_hat    <- IG_mlefunc(obs = y, w = covariates)
score.matrix <- IG_scorefunc(obs = y, mle = theta_hat, w = covariates)
pit.values   <- IG_pitfunc(obs = y , mle = theta_hat)

# Apply the goodness-of-fit test.
testYourModel(x = y, pit = pit.values, score = score.matrix)
```

    ## $Statistic
    ## Cramer-von-Mises Statistic 
    ##                  0.1402415 
    ## 
    ## $pvalue
    ## [1] 0.05200096

### References

\[1\] Stephens, M.A. (1974). \[EDF Statistics for Goodness of Fit and
Some Comparisons.\] *Journal of the American Statistical Association*,
Vol. 69, 730-737.

\[2\] Stephens, M.A. (1976). \[Asymptotic results for goodness-of-fit
statistics with unknown parameters.\] *Annals of Statistics*, Vol. 4,
357-369.
