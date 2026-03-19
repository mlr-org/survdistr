# Extract time points from a probability matrix or vector

Helper function to consistently obtain and validate the time points
associated with a probability matrix or vector.

## Usage

``` r
extract_times(x, times = NULL)
```

## Arguments

- x:

  ([`numeric()`](https://rdrr.io/r/base/numeric.html) \|
  [`matrix()`](https://rdrr.io/r/base/matrix.html))  
  Probability vector (length = time points) or matrix (rows =
  observations, columns = time points).

- times:

  ([`numeric()`](https://rdrr.io/r/base/numeric.html) \| `NULL`)  
  Optional vector of time points corresponding to `x`.

## Value

A validated numeric vector of time points.
