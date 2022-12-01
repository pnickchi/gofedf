README
================
2022-11-15

# GoFTest

<!-- badges: start -->

[![CRANstatus](https://www.r-pkg.org/badges/version/GoFTest)](https://cran.r-project.org/package=GoFTest)
[![GitHubRelease](https://img.shields.io/github/release/GoFTest/GoFTest?style=flat)](https://github.com/GoFTest/PEIMAN2/releases)
[![Github All Releases](https://img.shields.io/github/downloads/GoFTest/PEIMAN2/total.svg?style=flat)](https://github.com/GoFTest/PEIMAN2)

<!-- badges: end -->

Authors:

- [Payman Nickchi](https://github.com/pnickchi), <pnickchi@sfu.ca>
  (Maintainer)
- [Richard Lockhart](http://www.sfu.ca/~lockhart/), <lockhart@sfu.ca>

The GoFTest package provides tools to apply goodness of fit tests based
on empirical distribution function theory. The software provides
functions and routines to test the hypothesis that a sample follows a
distribution by calculating Cramer-von Mises or Anderson-Darling
statistic and computing the approximate pvalue by Imhof method.

The package can be used to apply goodness of fit test in any of the
following cases:

1)  Check if the residuals of a generalized linear model (with any link)
    follows the Gamma distribution.
2)  Check if the residuals of a linear model follows Normal
    distribution.
3)  Check if a bivariate continuous set of data follows any specific
    distribution that is defined by user.
4)  Apply goodness of fit test for bivariate Normal, Gamma, and
    Exponential distributions.

## Installation

You can install the released version of GoFTest from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages('GoFTest')
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github('pnickchi/GoFTest')
```

## Example

### Bivariate Normal distribution

In this example, we show how to apply goodness of fit test over a vector
of data and check if the data follows a normal distribution.

``` r
# Reproducible example
set.seed(1)

# Randomly generate some data
n <- 50
sim_data <- rnorm(n)

# Test if the data follows a normal distribution, calculate Cramer-von Mises statistic and approximate pvalue
testNormal(x = sim_data, method = 'cvm')
```

    ## $Statistic
    ## [1] 0.06711927
    ## 
    ## $pvalue
    ## [1] 0.2622589

``` r
# Test if the data follows a normal distribution, calculate Anderson-Darling statistic and approximate pvalue
testNormal(x = sim_data, method = 'ad')
```

    ## $Statistic
    ## [1] 0.4788945
    ## 
    ## $pvalue
    ## [1] 0.1274992

### Bivariate Gamma distribution

In this example, we show how to apply goodness of fit test over a vector
of data and check if the data follows a Gamma distribution.

``` r
# Reproducible example
set.seed(2)

# Randomly generate some data
n <- 50
sim_data <- rgamma(n, shape = 3)

# Test if the data follows a Gamma distribution, calculate Cramer-von Mises statistic and approximate pvalue
testGamma(x = sim_data, method = 'cvm')
```

    ## $Statistic
    ## [1] 0.0916967
    ## 
    ## $pvalue
    ## [1] 0.1676796

### Linear model with normal residuals

``` r
# Reproducible example
set.seed(3)

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

    ## $Statistic
    ## [1] 0.07322979
    ## 
    ## $pvalue
    ## [1] 0.1542016

``` r
# Or alternatively just pass 'myfit' object directly instead of X and y:
# myfit <- lm(y ~ X - 1, x = TRUE, y = TRUE)
# testLMNormal(fit = myfit)
```

### Generalized linear model with Gamma residuals

``` r
# Reproducible example
set.seed(4)


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
    ## [1] 0.02131014
    ## 
    ## $pvalue
    ## [1] 0.9278594
    ## 
    ## $converged
    ## [1] TRUE

### User defined distributions

One of the most important features of the package is to allow users to
apply the goodness of fit test based on empirical distribution functions
for a user defined distribution. This is helpful when the functions and
procedures for a distribution is not included in any R packages. For
example, imagine you have a sample of size
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
such as
![x\_{1},x\_{2},...,x\_{n}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_%7B1%7D%2Cx_%7B2%7D%2C...%2Cx_%7Bn%7D "x_{1},x_{2},...,x_{n}")
from a model with
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")
unknown parameters such as
![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta").
The package requires 3 elements to apply GoF for your model:

1)  A function to calculate score values for each observation and
    returns a (n x p) matrix. Each rows calculates the score for each
    observation. Alternatively you might have the score values as a
    matrix which is also acceptable.

2)  A function to calculate probability inverse transformation of data.
    This is usually a function that returns a vector with n elements
    where each element is
    ![F^{-1}(x\_{i})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F%5E%7B-1%7D%28x_%7Bi%7D%29 "F^{-1}(x_{i})").
    You can alternatively pass a vector of values as well.

3)  Maximum likelihood estimate of
    ![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta").

``` r
set.seed(123)
n <- 50
x <- rexp(n)

# Case one: score matrix and pit values are provided.
mle.theta     <- expMLE(x)
score.matrix  <- expScore(x, theta = mle.theta)
pit           <- expFx(x, theta = mle.theta)
testYourModel(x, score = score.matrix, Fx = pit, mle = mle.theta)
```

    ## $Statistic
    ## [1] 0.08503541
    ## 
    ## $pvalue
    ## [1] 0.4519658

``` r
# Case two: Two functions are provided to calculate score and pit.
mle.theta     <- expMLE(x)
testYourModel(x, score = expScore, Fx = expFx, mle = mle.theta)
```

    ## $Statistic
    ## [1] 0.08503541
    ## 
    ## $pvalue
    ## [1] 0.4519658
