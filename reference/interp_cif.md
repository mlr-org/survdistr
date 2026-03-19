# Interpolate CIF matrix

Interpolates cumulative incidence (CIF) functions (corresponding to one
competing event only) using left-continuous constant interpolation.

## Usage

``` r
interp_cif(x, times = NULL, eval_times = NULL, add_times = TRUE, check = TRUE)
```

## Arguments

- x:

  ([`matrix()`](https://rdrr.io/r/base/matrix.html))  
  CIF matrix (rows = observations, columns = time points).

- times:

  ([`numeric()`](https://rdrr.io/r/base/numeric.html) \| `NULL`)  
  Anchor time points. If `NULL`, extracted from `colnames(x)`.

- eval_times:

  ([`numeric()`](https://rdrr.io/r/base/numeric.html) \| `NULL`)  
  New time points at which to interpolate. Values need to be unique,
  increasing non-negative numbers. If `NULL`, the anchor `times` are
  used.

- add_times:

  (`logical(1)`)  
  If `TRUE`, attach `eval_times` as colnames in the output matrix.

- check:

  (`logical(1)`)  
  If `TRUE`, run input matrix validation checks using
  [`assert_prob()`](https://survdistr.mlr-org.com/reference/assert_prob.md);
  set to `FALSE` to skip checks (NOT recommended for external use).

## Value

Interpolated CIF matrix.
