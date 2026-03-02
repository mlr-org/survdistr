#' @title Assert probability matrix
#'
#' @description
 #' Validates that the input is a proper probability matrix representing either
#' a survival function (`"surv"`), cumulative distribution function (`"cdf"`),
#' or cumulative incidence function (`"cif"`).
#' Uses the internal Rcpp function `c_assert_prob_matrix()`.
#'
#' @details
#' The following conditions must hold:
#'
#' 1. The input `x` is a numeric matrix with no missing values.
#' 2. Time points (`times`) are numeric, non-negative, unique, and sorted.
#'    If not supplied, they are derived from `colnames(x)` (coerced to `numeric`).
#' 3. The number of time points equals the number of columns of `x`.
#' 4. All values are valid probabilities, i.e. lie in \eqn{[0,1]}.
#' 5. Each row is monotone:
#'    - `"surv"`: non-increasing survival curves, i.e. \eqn{S(t_i) \ge S(t_{i+1})}.
#'    - `"cdf"` / `"cif"`: non-decreasing functions, i.e. \eqn{F(t_i) \le F(t_{i+1})}.
#' 6. Boundary condition at `t = 0`:
#'    - `"surv"`: \eqn{S(0) = 1}.
#'    - `"cdf"` / `"cif"`: \eqn{F(0) = 0}.
#'
#' @param x (`matrix()`)\cr
#'  A probability matrix.
#'  Rows correspond to observations and columns correspond to time points.
#' @param times (`numeric()`|`NULL`)\cr
#'  Optional numeric vector of time points corresponding to the columns of `x`.
#'  If `numeric()`, these will be returned after the checks are performed.
#'  If `NULL` (default), the time points are extracted from `colnames(x)`.
#' @param type (`character(1)`)\cr
#'  Type of probability function: `"surv"` (default), `"cdf"`, or `"cif"`.
#'
#' @return
#' Invisibly returns the validated numeric vector of time points.
#' Throws an error if validation fails.
#'
#' @examples
#' x = matrix(data = c(1, 0.6, 0.4,
#'                     0.8, 0.8, 0.7),
#'            nrow = 2, ncol = 3, byrow = TRUE)
#'
#' # Explicitly provide time points
#' assert_prob_matrix(x, times = c(12, 34, 42), type = "surv")
#'
#' # Or use column names as time points
#' colnames(x) = c(12, 34, 42)
#' assert_prob_matrix(x)
#'
#' # check CDF
#' assert_prob_matrix(1 - x, type = "cdf")
#'
#' @export
assert_prob_matrix = function(x, times = NULL, type = "surv") {
  assert_matrix(x, mode = "numeric", any.missing = FALSE, min.rows = 1, min.cols = 1)
  type = assert_choice(type, c("surv", "cdf", "cif"))

  # extract and check times
  times = extract_times(x, times)

  # Check values and monotonicity via Rcpp code
  if (!c_assert_prob_matrix(x, type = type)) {
    msg = switch(type,
      "surv" = "Survival probabilities must be non-increasing and in [0,1].",
      "cdf"  = "CDF probabilities must be non-decreasing and in [0,1].",
      "cif"  = "CIF probabilities must be non-decreasing and in [0,1].")
    stop(msg)
  }

  # boundary at t = 0
  if (times[1] == 0) {
    # Do we need some tolerance here? i.e. 1 - .Machine$double.eps^0.9 still counts as 1?
    switch(type,
      "surv" = if (!all(x[, 1] == 1)) stop("At t = 0, survival S(0) must equal 1."),
      "cdf"  = if (!all(x[, 1] == 0)) stop("At t = 0, CDF(0) must equal 0."),
      "cif"  = if (!all(x[, 1] == 0)) stop("At t = 0, CIF(0) must equal 0.")
    )
  }

  invisible(times)
}

#' @title Assert probability vector
#' 
#' @description
#' Validates that the input is a proper probability vector representing either
#' a survival function (`"surv"`), cumulative distribution function (`"cdf"`),
#' or cumulative incidence function (`"cif"`).
#' 
#' @param x (`numeric()`)\cr
#'  Probability vector.
#' @param times (`numeric()`|`NULL`)\cr
#'  Optional time points corresponding to `x`. If `NULL`, extracted from `names(x)`.
#' @param type (`character(1)`)\cr
#'  Type of probability function: `"surv"` (default), `"cdf"`, or `"cif"`.
#' 
#' @return Invisibly returns validated numeric time points.
#' @export
assert_prob_vec = function(x, times = NULL, type = "surv") {
  assert_numeric(x, any.missing = FALSE, min.len = 1)
  type = assert_choice(type, c("surv", "cdf", "cif"))

  times = extract_times(x, times)

  if (any(x < 0 | x > 1)) {
    stop("Probabilities must lie in [0,1].")
  }

  # monotonicity
  diffs = diff(x)
  if (type == "surv" && any(diffs > 0)) {
    stop("Survival probabilities must be non-increasing.")
  }
  if (type %in% c("cdf", "cif") && any(diffs < 0)) {
    stop("CDF/CIF probabilities must be non-decreasing.")
  }

  # boundary at t = 0
  if (times[1] == 0) {
    if (type == "surv" && x[1] != 1) {
      stop("At t = 0, survival S(0) must equal 1.")
    }
    if (type %in% c("cdf", "cif") && x[1] != 0) {
      stop("At t = 0, CDF/CIF(0) must equal 0.")
    }
  }

  invisible(times)
}

#' @title Assert probability matrix or vector
#' 
#' @description 
#' Dispatches to `assert_prob_matrix()` or `assert_prob_vec()` based on the input type.
#' @param x (`numeric()` | `matrix()`)\cr
#'   Survival vector or matrix (rows = observations, columns = time points).
#' @param times (`numeric()` | `NULL`)\cr
#'   Original time points. If `NULL`, extracted from names/colnames.
#' @param type (`character(1)`)\cr
#'  Type of probability function: `"surv"` (default), `"cdf"`, or `"cif"`.
#' 
#' @return Invisibly returns validated numeric time points.
#' 
#' @noRd
assert_prob = function(x, times = NULL, type = "surv") {
  if (is.matrix(x)) {
    assert_prob_matrix(x, times, type = type)
  } else {
    assert_prob_vec(x, times, type = type)
  }
}

#' @title Assert Positive matrix
#'
#' @description Asserts if the given input is a non-negative matrix.
#' The following checks are performed:
#'
#' 1. All values are non-negative (i.e. \eqn{\ge 0}).
#' 2. Column names correspond to time points and should therefore be coercable to
#' `numeric` and increasing
#'
#' @param x (`matrix()`)\cr
#' Input matrix (rows = observations, columns = time points).
#'
#' @return if the assertion fails an error occurs, otherwise `NULL` is returned
#' invisibly.
#'
#' @examples
#' x = matrix(data = c(0,0.2,0.1,0.5,0.05,0.7), nrow = 2, ncol = 3, byrow = TRUE)
#' colnames(x) = c(12, 34, 42)
#' x
#'
#' assert_pos_mat(x)
#'
#' @noRd
assert_pos_mat = function(x) {
  assert_matrix(x, mode = "numeric", any.missing = FALSE, min.rows = 1, min.cols = 1, col.names = "named")
  assert_numeric(as.numeric(colnames(x)), unique = TRUE, sorted = TRUE,
                 any.missing = FALSE, null.ok = FALSE)
  stopifnot(x >= 0)
  invisible(NULL)
}
