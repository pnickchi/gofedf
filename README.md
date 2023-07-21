README
================
2023-07-21

# gofedf

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/gofedf)](https://cran.r-project.org/package=gofedf)
[![GitHub
Release](https://img.shields.io/github/release/pnickchi/gofedf?style=flat)](https://github.com/pnickchi/gofedf/releases)
[![Github All
Releases](https://img.shields.io/github/downloads/pnickchi/gofedf/total.svg?style=flat)](https://github.com/pnickchi/gofedf)
<!-- badges: end -->

Authors:

- [Richard Lockhart](http://www.sfu.ca/~lockhart/), <lockhart@sfu.ca>
- [Payman Nickchi](https://github.com/pnickchi), <pnickchi@sfu.ca>
  (Maintainer)

The `gofedf` package provides tools to apply goodness-of-fit tests based
on empirical distribution function theory. The package offers functions
and routines to test the hypothesis that a univariate sample follows a
distribution based on the empirical distribution function. The theory is
founded on reducing the problem to a stochastic process and computing
its covariance function. An approximate p-value, computed using the
farebrother method based on the limiting distribution of the statistic,
is provided. Users can run the test by calculating either the Cramer-von
Mises or Anderson-Darling statistic. The covariance function of the
stochastic process relies on specific characteristics of the assumed
model. Notably, knowledge of the Fisher information matrix and the
partial derivatives of the cumulative distribution function is crucial
for computing the covariance function. However, obtaining these
quantities can be computationally intensive or challenging in general
likelihood models. To overcome this limitation, we propose an
alternative method for estimating the covariance function of the
stochastic process directly from the sample data. The package provides
tools for this estimation and for testing if a sample comes from any
general likelihood model. In summary, the package can be used to apply a
goodness-of-fit test in any of the following settings:

- Verify if the assumptions about the response variable in a generalized
  linear model (with any link function) are satisfied. The current
  version only checks for the Gamma distribution.

- Validate if the residuals of a linear model follow a Normal
  distribution.

- Apply a goodness-of-fit test to examine if a set of bivariate samples
  follows a Normal, Gamma, or Exponential distribution.

- Determine if a bivariate continuous set of data conforms to a specific
  distribution defined by the user. In this case, probability inverted
  values of the sample and the score function (if any parameter
  estimation is involved) are required.

## Installation

You can install the released version of `gofedf` from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages('gofedf')
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github('pnickchi/gofedf')
```

## Examples

### 1. Bivariate Normal distribution

``` r
# Reproducible example
set.seed(123)

# Randomly generate some data
n <- 50
sim_data <- rnorm(n)

# Test if the data follows a Normal distribution, calculate Cramer-von Mises statistic and approximate p-value of the test.
testNormal(x = sim_data, method = 'cvm')
```

    ## $Statistic
    ## [1] 0.03781322
    ## 
    ## $pvalue
    ## [1] 0.6646717

``` r
# Test if the data follows a Normal distribution, calculate Anderson-Darling statistic and approximate p-value of the test.
testNormal(x = sim_data, method = 'ad')
```

    ## $Statistic
    ## [1] 0.2179704
    ## 
    ## $pvalue
    ## [1] 0.7833757

### 2. Bivariate Gamma distribution

``` r
# Reproducible example
set.seed(123)

# Randomly generate some data
n <- 50
sim_data <- rgamma(n, shape = 3)

# Test if the data follows a Gamma distribution, calculate Cramer-von Mises statistic and approximate pvalue
testGamma(x = sim_data, method = 'cvm')
```

    ## $Statistic
    ## [1] 0.0549759
    ## 
    ## $pvalue
    ## [1] 0.3318735

### 3. Linear model with normal residuals

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

# Compute response variable
y <- X %*% b + e

# Test if the residuals of the model follows a Normal distribution, calculate Cramer-von Mises statistic and approximate pvalue
testLMNormal(x = X, y)
```

    ## $Statistic
    ## [1] 0.02437418
    ## 
    ## $pvalue
    ## [1] 0.8740983

``` r
# Or alternatively just pass 'myfit' object directly instead of X and y:
# myfit <- lm(y ~ X - 1, x = TRUE, y = TRUE)
# testLMNormal(fit = myfit)
```

### 4. Generalized linear model with Gamma response and log link function

``` r
# Reproducible example
set.seed(123)


# Create a set of explanatory variables and response according to a generalized linear model with log link
n <- 50
p <- 2
X <- matrix( rnorm(n*p, mean = 10, sd = 0.1), nrow = n, ncol = p)
b <- runif(p)
e <- rgamma(n, shape = 3)
y <- exp(X %*% b) * e

# Test if the residuals of the model follows a Gamma distribution, calculate Cramer-von Mises statistic and approximate pvalue
testGLMGamma(x=X, y, l = 'log', method = 'cvm')
```

    ## $Statistic
    ## [1] 0.04026585
    ## 
    ## $pvalue
    ## [1] 0.5894211
    ## 
    ## $converged
    ## [1] TRUE

### 5. General likelihood model

One of the most important features of the package is to allow users to
apply the goodness-of-fit test based on empirical distribution functions
for any general likelihood model. This is helpful when the functions and
procedures for a distribution are not included in any R packages. For
example, consider you have a sample of size $n$, such as
$x_1, x_2, \ldots, x_n$, from a model with $p$ unknown parameters, such
as $\theta$. The package requires the followings to apply the test:

1)  A function to compute score function which returns a $n \times p$
    matrix. Alternatively you can pass the matrix directly.

2)  A vector of length $n$ contains the probability inverse transformed
    (PIT) values of sample, i.e $F^{-1}(x_{i})$.

In the following example, we show how to apply the test for Inverse
Gaussian distribution. Note that `IG_scorefunc` is a function to return
score function and `IG_pitfunc` returns a vector of PIT values.

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
