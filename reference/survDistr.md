# Survival Distribution Container

survDistr is an [R6](https://r6.r-lib.org/reference/R6Class.html)
specialized container designed for storing and managing prediction
outputs from survival models in single-event settings. This includes
models such as Cox proportional hazards, random survival forests, and
other classical or machine learning-based survival estimators.

The main prediction data type is survival matrix, where **rows represent
observations and columns represent time points**.

## Details

Key design features:

- The survival matrix is stored internally and can be accessed using the
  `$data()` method.

- The `$times` active field provides the time points corresponding to
  the matrix columns.

- The interpolation method is controlled via the `$method` active field.

- Survival-related quantities (e.g., distribution, density, hazard
  functions) are interpolated using the
  [`interp()`](https://survdistr.mlr-org.com/reference/interp.md)
  function.

- The
  [`assert_prob()`](https://survdistr.mlr-org.com/reference/assert_prob.md)
  function validates the input data matrix during construction if
  `check` is set to `TRUE`.

- Use the `$filter()` method to subset observations in-place by
  filtering rows of the stored matrix.

- Use `trim_duplicates = TRUE` in the constructor to remove flat
  survival segments (repeated values across time points) with a
  pre-specified tolerance (for a more controlled pre-processing, see
  [`trim_duplicates()`](https://survdistr.mlr-org.com/reference/trim_duplicates.md)).

## Active bindings

- `times`:

  (`numeric`)  
  Numeric vector of time points corresponding to columns of `data`.
  Read-only.

- `method`:

  (`character(1)`)  
  Interpolation method; one of `"const_surv"` (default), `"const_dens"`
  (alias: `"linear_surv"`) and `"const_haz"` (alias: `"exp_surv"`).

## Methods

### Public methods

- [`survDistr$new()`](#method-survDistr-new)

- [`survDistr$print()`](#method-survDistr-print)

- [`survDistr$data()`](#method-survDistr-data)

- [`survDistr$filter()`](#method-survDistr-filter)

- [`survDistr$survival()`](#method-survDistr-survival)

- [`survDistr$distribution()`](#method-survDistr-distribution)

- [`survDistr$density()`](#method-survDistr-density)

- [`survDistr$cumhazard()`](#method-survDistr-cumhazard)

- [`survDistr$hazard()`](#method-survDistr-hazard)

- [`survDistr$clone()`](#method-survDistr-clone)

------------------------------------------------------------------------

### Method `new()`

Creates a new instance of this
[R6](https://r6.r-lib.org/reference/R6Class.html) class.

#### Usage

    survDistr$new(
      x,
      times = NULL,
      method = "const_surv",
      check = TRUE,
      trim_duplicates = FALSE
    )

#### Arguments

- `x`:

  (`matrix`)  
  A numeric matrix of survival probabilities (values between 0 and 1).
  Column names must correspond to time points if `times` is `NULL`.

- `times`:

  (`numeric(1)`)  
  Numeric vector of time points for matrix `x`, must match the number of
  columns.

- `method`:

  (`character(1)`)  
  Interpolation method; one of `"const_surv"` (default), `"const_dens"`
  (alias: `"linear_surv"`) and `"const_haz"` (alias: `"exp_surv"`).

- `check`:

  (`logical(1)`)  
  If `TRUE`, run input matrix validation checks using
  [`assert_prob()`](https://survdistr.mlr-org.com/reference/assert_prob.md);
  set to `FALSE` to skip checks (NOT recommended for external use).

- `trim_duplicates`:

  (`logical(1)`)  
  If `TRUE`, removes adjacent duplicate values from the input using
  [`trim_duplicates()`](https://survdistr.mlr-org.com/reference/trim_duplicates.md).
  This eliminates flat segments in survival curves and improves
  interpolation efficiency. Default is `FALSE`.

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Displays summary information about a survDistr object, including the
number of observations and time points.

#### Usage

    survDistr$print()

------------------------------------------------------------------------

### Method [`data()`](https://rdrr.io/r/utils/data.html)

Return the stored data matrix.

#### Usage

    survDistr$data(rows = NULL, add_times = TRUE)

#### Arguments

- `rows`:

  ([`integer()`](https://rdrr.io/r/base/integer.html) \|
  [`logical()`](https://rdrr.io/r/base/logical.html) \| `NULL`)  
  Row indices or a logical vector used to filter observations. Logical
  vectors must have length equal to the number of observations. Integer
  indices must be positive and within range. If `NULL`, no filtering is
  applied.

- `add_times`:

  (`logical(1)`)  
  If `TRUE`, attach `eval_times` to the output.

#### Returns

(`matrix`)

------------------------------------------------------------------------

### Method [`filter()`](https://rdrr.io/r/stats/filter.html)

Filters observations **in-place** by subsetting rows of the stored
matrix.

#### Usage

    survDistr$filter(rows = NULL)

#### Arguments

- `rows`:

  ([`integer()`](https://rdrr.io/r/base/integer.html) \|
  [`logical()`](https://rdrr.io/r/base/logical.html) \| `NULL`)  
  Row indices or a logical vector used to filter observations. Logical
  vectors must have length equal to the number of observations. Integer
  indices must be positive and within range. If `NULL`, no filtering is
  applied.

#### Returns

Invisibly returns the `survDistr` object itself.

------------------------------------------------------------------------

### Method `survival()`

Computes survival probabilities \\S(t)\\ at the specified time points.

#### Usage

    survDistr$survival(rows = NULL, times = NULL, add_times = TRUE)

#### Arguments

- `rows`:

  ([`integer()`](https://rdrr.io/r/base/integer.html) \|
  [`logical()`](https://rdrr.io/r/base/logical.html) \| `NULL`)  
  Row indices or a logical vector used to filter observations. Logical
  vectors must have length equal to the number of observations. Integer
  indices must be positive and within range. If `NULL`, no filtering is
  applied.

- `times`:

  (`numeric`)  
  New time points at which to interpolate. Values need to be unique,
  increasing non-negative numbers. If `NULL`, the object's stored time
  points are used.

- `add_times`:

  (`logical(1)`)  
  If `TRUE`, attach `eval_times` to the output.

#### Returns

a `matrix` of survival probabilities (rows = observations, columns =
time points).

------------------------------------------------------------------------

### Method `distribution()`

Computes the cumulative distribution function \\F(t) = 1 - S(t)\\ or CDF
at the specified time points. \\F(t)\\ is the probability that the event
has occurred up until time \\t\\.

#### Usage

    survDistr$distribution(rows = NULL, times = NULL, add_times = TRUE)

#### Arguments

- `rows`:

  ([`integer()`](https://rdrr.io/r/base/integer.html) \|
  [`logical()`](https://rdrr.io/r/base/logical.html) \| `NULL`)  
  Row indices or a logical vector used to filter observations. Logical
  vectors must have length equal to the number of observations. Integer
  indices must be positive and within range. If `NULL`, no filtering is
  applied.

- `times`:

  (`numeric`)  
  New time points at which to interpolate. Values need to be unique,
  increasing non-negative numbers. If `NULL`, the object's stored time
  points are used.

- `add_times`:

  (`logical(1)`)  
  If `TRUE`, attach `eval_times` to the output.

#### Returns

a `matrix` of CDF values (rows = observations, columns = time points).

------------------------------------------------------------------------

### Method [`density()`](https://rdrr.io/r/stats/density.html)

Computes the probability density function \\f(t)\\ or PDF at the
specified time points. \\f(t) = \frac{d}{dt} F(t)\\ is the probability
of the event occurring at the specific time \\t\\.

#### Usage

    survDistr$density(rows = NULL, times = NULL, add_times = TRUE)

#### Arguments

- `rows`:

  ([`integer()`](https://rdrr.io/r/base/integer.html) \|
  [`logical()`](https://rdrr.io/r/base/logical.html) \| `NULL`)  
  Row indices or a logical vector used to filter observations. Logical
  vectors must have length equal to the number of observations. Integer
  indices must be positive and within range. If `NULL`, no filtering is
  applied.

- `times`:

  (`numeric`)  
  New time points at which to interpolate. Values need to be unique,
  increasing non-negative numbers. If `NULL`, the object's stored time
  points are used.

- `add_times`:

  (`logical(1)`)  
  If `TRUE`, attach `eval_times` to the output.

#### Returns

a `matrix` of PDF values (rows = observations, columns = time points).

------------------------------------------------------------------------

### Method `cumhazard()`

Computes the cumulative hazard (accumulated risk up to time \\t\\) at
the specified time points as \\H(t) = -log(S(t))\\.

#### Usage

    survDistr$cumhazard(rows = NULL, times = NULL, add_times = TRUE, eps = 1e-12)

#### Arguments

- `rows`:

  ([`integer()`](https://rdrr.io/r/base/integer.html) \|
  [`logical()`](https://rdrr.io/r/base/logical.html) \| `NULL`)  
  Row indices or a logical vector used to filter observations. Logical
  vectors must have length equal to the number of observations. Integer
  indices must be positive and within range. If `NULL`, no filtering is
  applied.

- `times`:

  (`numeric`)  
  New time points at which to interpolate. Values need to be unique,
  increasing non-negative numbers. If `NULL`, the object's stored time
  points are used.

- `add_times`:

  (`logical(1)`)  
  If `TRUE`, attach `eval_times` to the output.

- `eps`:

  (`numeric(1)`)  
  Small positive value used to replace extremely low survival
  probabilities when computing cumulative hazard, preventing numerical
  instability in \\-\log S(t)\\ calculations.

#### Returns

a `matrix` of cumulative hazards (rows = observations, columns = time
points).

------------------------------------------------------------------------

### Method `hazard()`

Computes the hazard \\h(t) = \frac{f(t)}{S(t)}\\ at the specified time
points. Hazard is the conditional instantaneous event rate at time \\t\\
given survival up to time \\t\\.

#### Usage

    survDistr$hazard(rows = NULL, times = NULL, add_times = TRUE)

#### Arguments

- `rows`:

  ([`integer()`](https://rdrr.io/r/base/integer.html) \|
  [`logical()`](https://rdrr.io/r/base/logical.html) \| `NULL`)  
  Row indices or a logical vector used to filter observations. Logical
  vectors must have length equal to the number of observations. Integer
  indices must be positive and within range. If `NULL`, no filtering is
  applied.

- `times`:

  (`numeric`)  
  New time points at which to interpolate. Values need to be unique,
  increasing non-negative numbers. If `NULL`, the object's stored time
  points are used.

- `add_times`:

  (`logical(1)`)  
  If `TRUE`, attach `eval_times` to the output.

#### Returns

a `matrix` of hazard values (rows = observations, columns = time
points).

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    survDistr$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# generate survival matrix
mat = matrix(data = c(1,0.6,0.4,0.8,0.8,0.7), nrow = 2,
             ncol = 3, byrow = TRUE)
times = c(12, 34, 42)
x = survDistr$new(mat, times)
x
#> A [2 x 3] survival matrix
#> Number of observations: 2
#> Number of time points: 3
#> Interpolation method: Piecewise Constant Survival 

# stored survival matrix
x$data()
#>       12  34  42
#> [1,] 1.0 0.6 0.4
#> [2,] 0.8 0.8 0.7

# interpolation method
x$method
#> [1] "const_surv"

# time points
x$times
#> [1] 12 34 42

eval_times = c(10, 30, 42, 50)
# S(t) at given time points (constant interpolation)
x$survival(times = eval_times)
#>      10  30  42  50
#> [1,]  1 1.0 0.4 0.4
#> [2,]  1 0.8 0.7 0.7
# same but with linear interpolation
x$method = "linear_surv"
x$survival(times = eval_times)
#>             10        30  42  50
#> [1,] 1.0000000 0.6727273 0.4 0.2
#> [2,] 0.8333333 0.8000000 0.7 0.6

# Cumulative hazard
x$cumhazard(times = eval_times)
#>             10        30        42        50
#> [1,] 0.0000000 0.3964153 0.9162907 1.6094379
#> [2,] 0.1823216 0.2231436 0.3566749 0.5108256

# Density
x$density(times = eval_times)
#>              10         30     42     50
#> [1,] 0.00000000 0.01818182 0.0250 0.0250
#> [2,] 0.01666667 0.00000000 0.0125 0.0125

# Hazard
x$hazard(times = eval_times)
#>        10         30         42         50
#> [1,] 0.00 0.02702703 0.06250000 0.12500000
#> [2,] 0.02 0.00000000 0.01785714 0.02083333
```
