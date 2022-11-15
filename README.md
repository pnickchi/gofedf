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

- Authors:  
  Richard Lockhart, Payman Nickchi (<pnickchi@sfu.ca>)

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
