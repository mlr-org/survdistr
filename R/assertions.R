#' Assert probability matrix
#'
#' Validates that the input is a proper probability matrix representing either
#' a survival function, cumulative distribution function (CDF), cumulative
#' incidence function (CIF), discrete hazard, or discrete density.
#'
#' @param x (`matrix()`)\cr
#'  A probability matrix (rows = observations, columns = time points).
#' @param times (`numeric()`|`NULL`)\cr
#'  Optional time points corresponding to `x`. If `NULL`, extracted from `colnames(x)`.
#' @param type (`character(1)`)\cr
#'  Type of probability function: `"surv"` (default), `"cdf"`, `"cif"`, `"haz"`, or `"dens"`.
#'
#' @return Invisibly returns the validated numeric time points.
#' @noRd
#' @keywords internal
assert_prob_matrix = function(x, times = NULL, type = "surv") {
  assert_matrix(x, mode = "numeric", any.missing = FALSE, min.rows = 1, min.cols = 1)
  type = assert_choice(type, c("surv", "cdf", "cif", "haz", "dens"))
  times = extract_times(x, times)

  # Range and monotonicity checks via Rcpp code
  if (!c_assert_prob_matrix(x, type = type)) {
    msg = switch(type,
      "surv" = "Survival probabilities must be non-increasing and in [0,1].",
      "cdf"  = "CDF probabilities must be non-decreasing and in [0,1].",
      "cif"  = "CIF probabilities must be non-decreasing and in [0,1].",
      "haz"  = "Hazard probabilities must be in [0,1].",
      "dens" = "Density probabilities must be in [0,1]."
    )
    stop(msg)
  }

  # boundary at t = 0
  if (times[1] == 0) {
    # Do we need some tolerance here? i.e. 1 - .Machine$double.eps^0.9 still counts as 1?
    switch(type,
      "surv" = if (!all(x[, 1] == 1)) stop("At t = 0, survival must equal 1."),
      "cdf"  = if (!all(x[, 1] == 0)) stop("At t = 0, CDF must equal 0."),
      "cif"  = if (!all(x[, 1] == 0)) stop("At t = 0, CIF must equal 0."),
      "haz"  = stop("For discrete hazards, times must be strictly positive."),
      "dens" = stop("For discrete densities, times must be strictly positive.")
    )
  }

  invisible(times)
}

#' @title Assert probability matrix or vector
#'
#' @description
#' Validates that the input is a proper probability matrix or vector representing either
#' a survival function, cumulative distribution function (CDF), cumulative
#' incidence function (CIF), discrete hazard, or discrete density.
#'
#' @details
#' The following conditions must hold:
#'
#' 1. The input `x` is a numeric matrix with no missing values.
#' 2. Time points (`times`) are numeric, non-negative, unique, and increasing.
#'    If not supplied, they are derived from `(col)names(x)` (coerced to `numeric`).
#' 3. All values are valid probabilities, i.e. lie in \eqn{[0,1]}.
#' 4. Each row is monotone:
#'    - `"surv"`: non-increasing survival curves, i.e. \eqn{S(t_i) \ge S(t_{i+1})}.
#'    - `"cdf"` / `"cif"`: non-decreasing functions, i.e. \eqn{F(t_i) \le F(t_{i+1})}.
#'    - `"haz"` / `"dens"`: no monotonicity requirement.
#' 5. Boundary condition at `t = 0`:
#'    - `"surv"`: \eqn{S(0) = 1}.
#'    - `"cdf"` / `"cif"`: \eqn{F(0) = 0}.
#'    - `"haz"` / `"dens"`: \eqn{t_1 > 0} (otherwise, nonzero hazard/density at `t = 0`
#'  implies \eqn{S(0) \neq 1})
#'
#' @param x (`numeric()` | `matrix()`)\cr
#'   Survival vector or matrix (rows = observations, columns = time points).
#' @param times (`numeric()` | `NULL`)\cr
#'   Original time points. If `NULL`, extracted from names/colnames.
#' @param type (`character(1)`)\cr
#'  Type of probability function: `"surv"` (default), `"cdf"`, `"cif"`, `"haz"`, or `"dens"`.
#'
#' @return Invisibly returns the validated numeric time points.
#' @examples
#' x = matrix(data = c(1, 0.6, 0.4,
#'                     0.8, 0.8, 0.7),
#'            nrow = 2, ncol = 3, byrow = TRUE)
#'
#' # Explicitly provide time points
#' assert_prob(x, times = c(12, 34, 42), type = "surv")
#'
#' # Or use column names as time points
#' colnames(x) = c(12, 34, 42)
#' assert_prob(x)
#'
#' # check CDF
#' assert_prob(1 - x, type = "cdf")
#'
#' # check discrete hazards
#' assert_prob(c(0.2, 0.01, 0.3), times = c(1, 2, 3), type = "haz")
#'
#' # check discrete densities
#' assert_prob(c(0.2, 0.01, 0.3), times = c(1, 2, 3), type = "dens")
#'
#' @export
assert_prob = function(x, times = NULL, type = "surv") {
  if (is.matrix(x)) {
    assert_prob_matrix(x, times, type = type)
  } else if (is.numeric(x)) {
    # vector input: coerce to matrix with colnames if possible
    mat = matrix(x, nrow = 1)
    vec_names = names(x)
    if (!is.null(vec_names)) colnames(mat) = vec_names

    assert_prob_matrix(mat, times, type = type)
  } else {
    stop("Must be of type numeric or matrix.")
  }
}
