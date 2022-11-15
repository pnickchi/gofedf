README
================
2022-11-14

# GoFTest

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/PEIMAN2)](https://cran.r-project.org/package=PEIMAN2)
[![GitHub
Release](https://img.shields.io/github/release/jafarilab/PEIMAN2?style=flat)](https://github.com/jafarilab/PEIMAN2/releases)
[![Github All
Releases](https://img.shields.io/github/downloads/jafarilab/PEIMAN2/total.svg?style=flat)](https://github.com/jafarilab/PEIMAN2)

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
