# Assert probability matrix or vector

Validates that the input is a proper probability matrix or vector
representing either a survival function, cumulative distribution
function (CDF), cumulative incidence function (CIF), discrete hazard, or
discrete density.

## Usage

``` r
assert_prob(x, times = NULL, type = "surv")
```

## Arguments

- x:

  ([`numeric()`](https://rdrr.io/r/base/numeric.html) \|
  [`matrix()`](https://rdrr.io/r/base/matrix.html))  
  Survival vector or matrix (rows = observations, columns = time
  points).

- times:

  ([`numeric()`](https://rdrr.io/r/base/numeric.html) \| `NULL`)  
  Original time points. If `NULL`, extracted from names/colnames.

- type:

  (`character(1)`)  
  Type of probability function: `"surv"` (default), `"cdf"`, `"cif"`,
  `"haz"`, or `"dens"`.

## Value

Invisibly returns the validated numeric time points.

## Details

The following conditions must hold:

1.  The input `x` is a numeric matrix with no missing values.

2.  Time points (`times`) are numeric, non-negative, unique, and
    increasing. If not supplied, they are derived from `(col)names(x)`
    (coerced to `numeric`).

3.  All values are valid probabilities, i.e. lie in \\\[0,1\]\\.

4.  Each row is monotone:

    - `"surv"`: non-increasing survival curves, i.e. \\S(t_i) \ge
      S(t\_{i+1})\\.

    - `"cdf"` / `"cif"`: non-decreasing functions, i.e. \\F(t_i) \le
      F(t\_{i+1})\\.

    - `"haz"` / `"dens"`: no monotonicity requirement.

5.  Boundary condition at `t = 0`:

    - `"surv"`: \\S(0) = 1\\.

    - `"cdf"` / `"cif"`: \\F(0) = 0\\.

    - `"haz"` / `"dens"`: \\t_1 \> 0\\ (otherwise, nonzero
      hazard/density at `t = 0` implies \\S(0) \neq 1\\)

## Examples

``` r
x = matrix(data = c(1, 0.6, 0.4,
                    0.8, 0.8, 0.7),
           nrow = 2, ncol = 3, byrow = TRUE)

# Explicitly provide time points
assert_prob(x, times = c(12, 34, 42), type = "surv")

# Or use column names as time points
colnames(x) = c(12, 34, 42)
assert_prob(x)

# check CDF
assert_prob(1 - x, type = "cdf")

# check discrete hazards
assert_prob(c(0.2, 0.01, 0.3), times = c(1, 2, 3), type = "haz")

# check discrete densities
assert_prob(c(0.2, 0.01, 0.3), times = c(1, 2, 3), type = "dens")
```
