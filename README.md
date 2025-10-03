
# survdistr

⚠️ **Under Development**  
This package is experimental and not yet intended for general use.  
APIs may change without notice, and functionality is incomplete.  
Please use only for testing and development purposes.

Survival distribution containers for efficient storage, management, and
evaluation of survival model predictions.

<!-- badges: start -->

[![r-cmd-check](https://github.com/mlr-org/survdistr/actions/workflows/r-cmd-check.yml/badge.svg)](https://github.com/mlr-org/mlr3survival/actions/workflows/r-cmd-check.yml)
[![Codecov test
coverage](https://codecov.io/gh/mlr-org/survdistr/branch/main/graph/badge.svg)](https://app.codecov.io/gh/mlr-org/survdistr?branch=main)
<!--[![CRAN Status](https://www.r-pkg.org/badges/version-ago/survdistr)](https://cran.r-project.org/package=survdistr) -->
<!-- badges: end -->

## Installation

Install the development version from GitHub:

``` r
# install.packages("pak")
pak::pak("mlr-org/survdistr")
```

## Example

Linear interpolation of a survival matrix using the `survDistr` R6
class:

``` r
library(survdistr)

# generate survival matrix
mat = matrix(data = c(1,0.6,0.4,0.8,0.8,0.7), nrow = 2,
             ncol = 3, byrow = TRUE)
x = survDistr$new(x = mat, times = c(12, 34, 42), interp_meth = "linear_surv")
x
```

    ## A [2 x 3] survival matrix
    ## Number of observations: 2
    ## Number of time points: 3
    ## Interpolation method: Piece-wise Linear Survival

``` r
# stored survival matrix
x$data()
```

    ##       12  34  42
    ## [1,] 1.0 0.6 0.4
    ## [2,] 0.8 0.8 0.7

``` r
# S(t) at given time points (linear interpolation)
x$survival(times = c(5, 30, 42, 50))
```

    ##              5        30  42  50
    ## [1,] 1.0000000 0.6727273 0.4 0.2
    ## [2,] 0.9166667 0.8000000 0.7 0.6

``` r
# Cumulative hazard H(t)
x$cumhazard(times = c(50, 5, 5, 42)) # times can be unordered or duplicated
```

    ##             50          5          5        42
    ## [1,] 1.6094379 0.00000000 0.00000000 0.9162907
    ## [2,] 0.5108256 0.08701138 0.08701138 0.3566749

Interpolation of a Kaplan-Meier survival curve using exported R function
that calls C code:

``` r
library(survival)
library(data.table)

fit = survfit(formula = Surv(time, status) ~ 1, data = veteran)
tab = data.table(time = fit$time, surv = fit$surv)
tab
```

    ##       time        surv
    ##      <num>       <num>
    ##   1:     1 0.985401460
    ##   2:     2 0.978102190
    ##   3:     3 0.970802920
    ##   4:     4 0.963503650
    ##   5:     7 0.941605839
    ##  ---                  
    ##  97:   467 0.036018043
    ##  98:   553 0.027013532
    ##  99:   587 0.018009021
    ## 100:   991 0.009004511
    ## 101:   999 0.000000000

``` r
# constant S(t) interpolation
vec_interp(
  x = tab$surv, 
  times = tab$time, 
  eval_times = c(0, 3.5, 995, 1000),
  constant = TRUE,
  type = "surv"
)
```

    ##           0         3.5         995        1000 
    ## 1.000000000 0.970802920 0.009004511 0.000000000

``` r
# linear S(t) interpolation
vec_interp(
  x = tab$surv, 
  times = tab$time, 
  eval_times = c(0, 3.5, 995, 1000),
  constant = FALSE,
  type = "surv"
)
```

    ##           0         3.5         995        1000 
    ## 1.000000000 0.967153285 0.004502255 0.000000000
