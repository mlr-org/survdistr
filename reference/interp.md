# Interpolate Survival Curves

Interpolates survival curves (vector or matrix) at new time points using
internal C++ interpolation functions. Output can be survival, cumulative
distribution, density, hazard or cumulative hazard functions.

## Usage

``` r
interp(
  x,
  times = NULL,
  eval_times = NULL,
  method = "const_surv",
  output = "surv",
  add_times = TRUE,
  check = TRUE,
  eps = 1e-12,
  trim_duplicates = FALSE
)
```

## Arguments

- x:

  ([`numeric()`](https://rdrr.io/r/base/numeric.html) \|
  [`matrix()`](https://rdrr.io/r/base/matrix.html))  
  Survival vector or matrix (rows = observations, columns = time
  points).

- times:

  ([`numeric()`](https://rdrr.io/r/base/numeric.html) \| `NULL`)  
  Anchor time points. If `NULL`, extracted from names/colnames.

- eval_times:

  ([`numeric()`](https://rdrr.io/r/base/numeric.html) \| `NULL`)  
  New time points at which to interpolate. Values need to be unique,
  increasing non-negative numbers. If `NULL`, the anchor `times` are
  used.

- method:

  (`character(1)`)  
  Interpolation method; one of `"const_surv"` (default), `"const_dens"`
  (alias: `"linear_surv"`) and `"const_haz"` (alias: `"exp_surv"`).

- output:

  (`character(1)`)  
  Output type: `"surv"`, `"cdf"`, `"cumhaz"`, `"density"` or `"hazard"`.

- add_times:

  (`logical(1)`)  
  If `TRUE`, attach `eval_times` as names/colnames.

- check:

  (`logical(1)`)  
  If `TRUE`, run input matrix validation checks using
  [`assert_prob()`](https://survdistr.mlr-org.com/reference/assert_prob.md);
  set to `FALSE` to skip checks (NOT recommended for external use).

- eps:

  (`numeric(1)`)  
  Small positive value used to replace extremely low survival
  probabilities when computing cumulative hazard, preventing numerical
  instability in \\-\log S(t)\\ calculations.

- trim_duplicates:

  (`logical(1)`)  
  If `TRUE`, removes adjacent duplicate values from the input using
  [`trim_duplicates()`](https://survdistr.mlr-org.com/reference/trim_duplicates.md).
  This eliminates flat segments in survival curves and improves
  interpolation efficiency. Default is `FALSE`.

## Value

A numeric vector or matrix of interpolated values.

## Details

Input must always be *survival probabilities*. We currently provide
three interpolation options:

- `"const_surv"`: left-continuous constant interpolation of S(t)
  (default).

- `"const_dens"`/`"linear_surv"`: linear interpolation of S(t)
  (equivalent to piecewise constant interpolation of the density
  function).

- `"const_haz"`/`"exp_surv"`: exponential interpolation of S(t)
  (equivalent to piecewise constant interpolation of the hazard
  function).

For formulas for each method, see respective Tables in arxiv preprint
(TODO: add link).

For constant hazard interpolation (`"const_haz"`), any right-anchor S(t)
values equal to 0 are internally floored at `min(1e-12, S_left)` within
each interval. This keeps hazards/densities finite without allowing a
local increase in S(t).

## Examples

``` r
x = matrix(c(1, 0.8, 0.6,
             1, 0.7, 0.4), nrow = 2, byrow = TRUE)
times = c(0, 8, 13)
eval_times = c(5, 10, 14)

# constant S(t) interpolation
interp(x, times, eval_times)
#>      5  10  14
#> [1,] 1 0.8 0.6
#> [2,] 1 0.7 0.4

# linear S(t) interpolation
interp(x, times, eval_times, method = "linear_surv")
#>           5   10   14
#> [1,] 0.8750 0.72 0.56
#> [2,] 0.8125 0.58 0.34

# exponential S(t) interpolation (same as `method = "const_haz"`)
interp(x, times, eval_times, method = "exp_surv")
#>              5        10        14
#> [1,] 0.8698237 0.7130410 0.5664525
#> [2,] 0.8001774 0.5596066 0.3576452

# Cumulative distribution with linear S(t) interpolation
interp(x, times, eval_times, method = "linear_surv", output = "cdf")
#>           5   10   14
#> [1,] 0.1250 0.28 0.44
#> [2,] 0.1875 0.42 0.66

# H(t) with linear S(t) interpolation
interp(x, times, eval_times, method = "linear_surv", output = "cumhaz")
#>              5        10        14
#> [1,] 0.1335314 0.3285041 0.5798185
#> [2,] 0.2076394 0.5447272 1.0788097

# f(t) with constant hazard interpolation
interp(x, times, eval_times, method = "const_haz", output = "density")
#>               5         10         14
#> [1,] 0.02426194 0.04102582 0.03259165
#> [2,] 0.03567540 0.06263294 0.04002878

# h(t) with constant hazard interpolation
interp(x, times, eval_times, method = "const_haz", output = "hazard")
#>               5         10         14
#> [1,] 0.02789294 0.05753641 0.05753641
#> [2,] 0.04458437 0.11192316 0.11192316
```
