README
================
2023-07-18

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

In this example, we demonstrate how to apply a goodness-of-fit test
based on the empirical distribution function (EDF) to assess whether a
vector of data follows a Normal distribution.

``` r
# Reproducible example
set.seed(123)

# Randomly generate some data
n <- 50
sim_data <- rnorm(n)

# Test if the data follows a Normal distribution, calculate Cramer-von Mises statistic and approximate p-value of the test.
testNormal(x = sim_data, method = 'cvm')
```

    ## Warning in getpvalue(u = cvm, eigen = ev): CompQuadForm failed to generate a
    ## correct pvalue. The pvalue lies between 0.26737599937446 and 1

    ## $Statistic
    ## [1] 0.03781322
    ## 
    ## $pvalue
    ## [1] 8

``` r
# Test if the data follows a Normal distribution, calculate Anderson-Darling statistic and approximate p-value of the test.
testNormal(x = sim_data, method = 'ad')
```

    ## Warning in getpvalue(u = AD, eigen = ev): CompQuadForm failed to generate a
    ## correct pvalue. The pvalue lies between 0.263880190202964 and 1

    ## $Statistic
    ## [1] 0.2179704
    ## 
    ## $pvalue
    ## [1] 8

### 2. Bivariate Gamma distribution

In this example, we show how to apply goodness of fit test over a vector
of data and check if the data follows a Gamma distribution.

``` r
# Reproducible example
set.seed(123)

# Randomly generate some data
n <- 50
sim_data <- rgamma(n, shape = 3)

# Test if the data follows a Gamma distribution, calculate Cramer-von Mises statistic and approximate pvalue
testGamma(x = sim_data, method = 'cvm')
```

    ## Warning in getpvalue(u = cvm, eigen = ev): CompQuadForm failed to generate a
    ## correct pvalue. The pvalue lies between 0.144766470183369 and 0.993008632555107

    ## $Statistic
    ## [1] 0.0549759
    ## 
    ## $pvalue
    ## [1] 8

### 3. Linear model with normal residuals

``` r
# Reproducible example
set.seed(123)

# Create a set of explanatory variables and response according to a linear model
n <- 50
p <- 5
X <- matrix( runif(n*p), nrow = n, ncol = p)
e <- runif(n)
b <- runif(p)
y <- X %*% b + e

# Test if the residuals of the model follows a Normal distribution, calculate Cramer-von Mises statistic and approximate pvalue
testLMNormal(x = X, y)
```

    ## Warning in getpvalue(u = cvm, eigen = ev): CompQuadForm failed to generate a
    ## correct pvalue. The pvalue lies between 0.0144740066675692 and
    ## 0.341840751142698

    ## $Statistic
    ## [1] 0.0971424
    ## 
    ## $pvalue
    ## [1] 8

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

    ## Warning in getpvalue(u = cvm, eigen = ev): CompQuadForm failed to generate a
    ## correct pvalue. The pvalue lies between 0.258929317433151 and 1

    ## $Statistic
    ## [1] 0.04026585
    ## 
    ## $pvalue
    ## [1] 8
    ## 
    ## $converged
    ## [1] TRUE

### 4. General likelihood model

One of the most important features of the package is to allow users to
apply the goodness of fit test based on empirical distribution functions
for any general likelihood model. This is helpful when the functions and
procedures for a distribution is not included in any R packages. For
example, imagine you have a sample of size $n$ such as
$x_{1},x_{2},...,x_{n}$ from a model with $p$ unknown parameters such as
$\theta$. The package requires two inputs to apply the test:

1)  A function to calculate score values for each observation and
    returns a (n x p) matrix. Each rows calculates the score for each
    observation. Alternatively you might have the score values as a
    matrix which is also acceptable.

2)  A vector of length $n$ where elements where each element is
    $F^{-1}(x_{i})$.

``` r
# Example: Inverse Gaussian (IG) distribution with weights

# Reproducible example
#set.seed(123)


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

    ## Warning in getpvalue(u = cvm, eigen = ev): CompQuadForm failed to generate a
    ## correct pvalue. The pvalue lies between 0.00530976850160436 and
    ## 0.164692390389004

    ## $Statistic
    ## Cramer-von-Mises Statistic 
    ##                  0.1698982 
    ## 
    ## $pvalue
    ## [1] 8

### References

\[1\] Stephens, M.A. (1974). \[EDF Statistics for Goodness of Fit and
Some Comparisons.\] *Journal of the American Statistical Association*,
Vol. 69, 730-737. \[2\] Stephens, M.A. (1976). \[Asymptotic results for
goodness-of-fit statistics with unknown parameters.\] *Annals of
Statistics*, Vol. 4, 357-369.
