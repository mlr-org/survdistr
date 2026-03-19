
# survdistr

Survival distribution container for efficient storage, management, and
interpolation of survival model predictions.

<!-- badges: start -->

[![r-cmd-check](https://github.com/mlr-org/survdistr/actions/workflows/r-cmd-check.yml/badge.svg)](https://github.com/mlr-org/mlr3survival/actions/workflows/r-cmd-check.yml)
[![Codecov test
coverage](https://codecov.io/gh/mlr-org/survdistr/branch/main/graph/badge.svg)](https://app.codecov.io/gh/mlr-org/survdistr?branch=main)
[![CRAN
Status](https://www.r-pkg.org/badges/version-ago/survdistr)](https://cran.r-project.org/package=survdistr)
<!-- badges: end -->

## Installation

Install the development version from GitHub:

``` r
# install.packages("pak")
pak::pak("mlr-org/survdistr")
```

## Examples

Linear interpolation of a survival matrix using the `survDistr` R6
class:

``` r
library(survdistr)

# generate survival matrix
mat = matrix(data = c(0.9,0.6,0.4,0.8,0.8,0.7), nrow = 2,
             ncol = 3, byrow = TRUE)
x = survDistr$new(x = mat, times = c(12, 34, 42), method = "linear_surv")
x
```

    ## A [2 x 3] survival matrix
    ## Number of observations: 2
    ## Number of time points: 3
    ## Interpolation method: Piecewise Constant Density (Linear Survival)

``` r
# stored survival matrix
x$data()
```

    ##       12  34  42
    ## [1,] 0.9 0.6 0.4
    ## [2,] 0.8 0.8 0.7

``` r
# S(t) at requested time points (linear interpolation)
x$survival(times = c(5, 30, 42, 50))
```

    ##              5        30  42  50
    ## [1,] 0.9583333 0.6545455 0.4 0.2
    ## [2,] 0.9166667 0.8000000 0.7 0.6

``` r
# Cumulative hazard H(t)
x$cumhazard(times = c(5, 42))
```

    ##               5        42
    ## [1,] 0.04255961 0.9162907
    ## [2,] 0.08701138 0.3566749

``` r
# Probability density f(t)
x$density(times = c(5, 30, 42))
```

    ##                5         30     42
    ## [1,] 0.008333333 0.01363636 0.0250
    ## [2,] 0.016666667 0.00000000 0.0125

``` r
# Hazard h(t)
x$hazard(times = c(5, 30, 42))
```

    ##                5         30         42
    ## [1,] 0.008695652 0.02083333 0.06250000
    ## [2,] 0.018181818 0.00000000 0.01785714

Interpolation of a Kaplan-Meier survival curve using exported R function
that calls C++ code:

``` r
library(survival)

fit = survfit(formula = Surv(time, status) ~ 1, data = veteran)
tab = data.frame(time = fit$time, surv = fit$surv)
head(tab)
```

    ##   time      surv
    ## 1    1 0.9854015
    ## 2    2 0.9781022
    ## 3    3 0.9708029
    ## 4    4 0.9635036
    ## 5    7 0.9416058
    ## 6    8 0.9124088

``` r
tail(tab)
```

    ##     time        surv
    ## 96   411 0.045022553
    ## 97   467 0.036018043
    ## 98   553 0.027013532
    ## 99   587 0.018009021
    ## 100  991 0.009004511
    ## 101  999 0.000000000

``` r
# constant S(t) interpolation
interp(
  x = tab$surv,
  times = tab$time,
  eval_times = c(0, 3.5, 995)
)
```

    ##           0         3.5         995 
    ## 1.000000000 0.970802920 0.009004511

``` r
# linear S(t) interpolation
interp(
  x = tab$surv,
  times = tab$time,
  eval_times = c(0, 3.5, 995),
  method = "linear_surv"
)
```

    ##           0         3.5         995 
    ## 1.000000000 0.967153285 0.004502255

``` r
# exponential S(t) interpolation
interp(
  x = tab$surv,
  times = tab$time,
  eval_times = c(0, 3.5, 995),
  method = "exp_surv"
)
```

    ##         0       3.5       995 
    ## 1.0000000 0.9671464 0.0000000

## Code of Conduct

Please note that the survdistr project is released with a [Contributor
Code of Conduct](https://survdistr.mlr-org.com/CODE_OF_CONDUCT.html). By
contributing to this project, you agree to abide by its terms.
