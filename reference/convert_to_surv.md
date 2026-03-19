# Convert density/hazard to survival

Converts density or hazards from one of four input representations to
survival probabilities at the same anchor time points (no
interpolation).

## Usage

``` r
convert_to_surv(
  x,
  times = NULL,
  input = "cont_haz",
  check = TRUE,
  clamp_surv = FALSE,
  eps = 1e-12
)
```

## Arguments

- x:

  ([`numeric()`](https://rdrr.io/r/base/numeric.html) \|
  [`matrix()`](https://rdrr.io/r/base/matrix.html))  
  Input vector or matrix (rows = observations, columns = time points).

- times:

  ([`numeric()`](https://rdrr.io/r/base/numeric.html) \| `NULL`)  
  Anchor time points. If `NULL`, extracted from names/colnames of `x`.

- input:

  (`character(1)`)  
  Input type. One of `"disc_haz"`, `"disc_dens"`, `"cont_haz"` or
  `"cont_dens"`.

- check:

  (`logical(1)`)  
  If `TRUE` (default), run *input* validation checks. Disable only if
  you know the input is valid and want to skip checks for speed.

- clamp_surv:

  (`logical(1)`)  
  If `TRUE`, clamp survival probabilities to `[eps, 1]` to avoid
  numerical issues.

- eps:

  (`numeric(1)`)  
  Small value used to clamp near-zero survival probabilities if
  `clamp_surv = TRUE`.

## Value

A numeric vector or matrix of survival probabilities with the same
dimensions as `x`.

## Details

Let \\t_1,\dots,t_B\\ denote the anchor time points, \\\Delta_j = t_j -
t\_{j-1}\\, and \\S_j = S(t_j)\\ the survival probabilities at the
anchors. The conversion depends on the value of `input`:

- `"disc_dens"`: \\S_j = 1 - \sum\_{k=1}^j \tilde f_k\\

- `"disc_haz"`: \\S_j = \prod\_{k=1}^j (1 - \tilde h_k)\\

- `"cont_dens"`: \\S_j = 1 - \sum\_{k=1}^j f_k \Delta_k\\

- `"cont_haz"`: \\S_j = \exp\\\left(-\sum\_{k=1}^j \lambda_k
  \Delta_k\right)\\

## Validation

If `check = TRUE`, we validate that the input is a proper discrete
density/hazard matrix or vector using
[`assert_prob()`](https://survdistr.mlr-org.com/reference/assert_prob.md).
For continuous hazards/densities, we only check that the input is a
non-negative numeric matrix/vector.

## Examples

``` r
# Continuous hazard => survival
haz_cont = c(0.02, 0.1, 0.2, 0.15)
times = c(0, 1, 2, 3)
convert_to_surv(haz_cont, times = times, input = "cont_haz")
#> [1] 1.0000000 0.9048374 0.7408182 0.6376282

# Discrete hazard => survival
haz_disc = c(0.1, 0.2, 0.15)
times = c(1, 2, 3)
convert_to_surv(haz_disc, times = times, input = "disc_haz")
#> [1] 0.900 0.720 0.612
```
