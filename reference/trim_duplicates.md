# Remove adjacent duplicate values

Removes adjacent duplicate values over the time axis, possibly from a
probability vector or matrix (e.g. survival curves). Equality is
determined with a numeric tolerance.

For matrices, duplicate detection is done column-wise across all rows.
Only the earliest time point in each run of (near-)equal values is kept.

## Usage

``` r
trim_duplicates(x, times = NULL, tol = 1e-10)
```

## Arguments

- x:

  ([`numeric()`](https://rdrr.io/r/base/numeric.html) \|
  [`matrix()`](https://rdrr.io/r/base/matrix.html)) Vector (length =
  time points) or matrix (rows = observations, columns = time points).

- times:

  ([`numeric()`](https://rdrr.io/r/base/numeric.html) \| `NULL`)
  Optional time points corresponding to `x`. If `NULL`, extracted from
  names/colnames.

- tol:

  (`numeric(1)`) Absolute tolerance used to detect equality between
  adjacent time points.

## Value

A named list with:

- `x`: numeric vector or matrix with duplicate adjacent time points
  removed.

- `times`: numeric vector of retained time points.

## Examples

``` r
# remove adjacent duplicates from a survival probability vector
surv = c(1, 1, 0.8, 0.8, 0.5, 0.5, 0.2)
trim_duplicates(surv, times = 1:7)
#> $x
#> [1] 1.0 0.8 0.5 0.2
#> 
#> $times
#> [1] 1 3 5 7
#> 
```
