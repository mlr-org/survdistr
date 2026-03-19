# Coerce Object to [survDistr](https://survdistr.mlr-org.com/reference/survDistr.md)

S3 generic to coerce supported objects into a
[survDistr](https://survdistr.mlr-org.com/reference/survDistr.md)
object.

## Usage

``` r
as_survDistr(x, ...)

# S3 method for class 'matrix'
as_survDistr(
  x,
  times = NULL,
  method = "const_surv",
  check = TRUE,
  trim_duplicates = FALSE,
  ...
)

# S3 method for class 'survDistr'
as_survDistr(x, ...)

# Default S3 method
as_survDistr(x, ...)
```

## Arguments

- x:

  Object to convert.

- ...:

  Additional arguments passed to methods.

- times:

  (`numeric`) Numeric vector of time points corresponding to columns of
  `x`. If `NULL`, column names of `x` are used.

- method:

  (`character(1)`)  
  Interpolation method passed to
  [survDistr](https://survdistr.mlr-org.com/reference/survDistr.md)
  constructor.

- check:

  (`logical(1)`)  
  Whether to validate `x` and `times`.

- trim_duplicates:

  (`logical(1)`)  
  Whether to remove duplicate S(t) values and corresponding time points.

## Value

A [survDistr](https://survdistr.mlr-org.com/reference/survDistr.md)
object.
