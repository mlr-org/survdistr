#' @title Assert probability matrix
#'
#' @description
#' Asserts whether the given input is a valid (discrete) survival/CDF/CIF probability
#' matrix using the internal Rcpp function `c_assert_prob_matrix()`.
#'
#' The following checks are performed:
#'
#' 1. The input `x` is a numeric matrix with no missing values.
#' 2. Time points (`times`) are numeric, non-negative, unique, and sorted.
#'    If not supplied, they are derived from the column names of `x`
#'    (coerced to numeric, possibly with loss of accuracy).
#' 3. The number of time points matches the number of columns of `x`.
#' 4. All values are valid probabilities, i.e. \eqn{S(t) \in [0,1]}.
#' 5. Per row (observation-wise), the survival probabilities decrease
#'    non-strictly, i.e. \eqn{S(t_i) \ge S(t_{i+1})}.
#' 6. If the first time point is 0, then \eqn{S(t = 0) = 1}.
#'
#' @param x (`matrix()`)\cr
#'  A probability matrix.
#'  Rows correspond to observations and columns correspond to time points.
#' @param times (`numeric()`|`NULL`)\cr
#'  Optional numeric vector of time points corresponding to the columns of `x`.
#'  If `numeric()`, these will be returned after the checks are performed.
#'  If `NULL`, the time points are extracted from `colnames(x)`.
#' @param type (`character(1)`)\cr
#'  One of `"surv"`, `"cdf"`, `"cif"`.
#'
#' @return
#' If the assertion fails, an error is thrown.
#' Otherwise, the validated numeric vector of time points is returned invisibly.
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
#' @export
assert_prob_matrix = function(x, times = NULL, type = "surv") {
  assert_matrix(x, any.missing = FALSE, min.rows = 1, min.cols = 1)
  type = assert_choice(type, c("surv", "cdf", "cif"))

  # check or derive times
  if (is.null(times)) {
    if (is.null(colnames(x))) {
      stop("Column names are required if 'times' is not provided.")
    }
    times = assert_numeric(as.numeric(colnames(x)),
                           lower = 0, unique = TRUE, sorted = TRUE,
                           any.missing = FALSE)
  } else {
    times = assert_numeric(times,
                           lower = 0, unique = TRUE, sorted = TRUE,
                           any.missing = FALSE, null.ok = FALSE)
    if (ncol(x) != length(times)) {
      stop("Number of columns in 'x' must match length of provided 'times'.")
    }
  }

  # Check values and monotonicity via Rcpp code
  if (!c_assert_prob_matrix(x, type = type)) {
    msg = switch(type,
                 "surv" = "Survival probabilities must be non-increasing and in [0,1].",
                 "cdf"  = "CDF probabilities must be non-decreasing and in [0,1].",
                 "cif"  = "CIF probabilities must be non-decreasing and in [0,1].")
    stop(msg)
  }

  if (times[1] == 0) {
    if (type == "surv") {
      if (!all(x[, 1] == 1)) stop("S(0) must equal 1.")
    } else {
      # type is either "cdf" or "cif"
      if (!all(x[, 1] == 0)) stop("CDF/CIF(0) must equal 0.")
    }
  }

  invisible(times)
}

#' @title Assert hazard matrix
#'
#' @description Asserts if the given input is a (discrete) hazard matrix.
#' The following checks are performed:
#'
#' 1. All values are non-negative, i.e. \eqn{h(t) \ge 0}
#' 2. Column names correspond to time points and should therefore be coercable to
#' `numeric` and increasing
#'
#' @param x (`matrix()`)\cr
#' A matrix of hazards.
#' Rows are observations and columns are (increasing) time points.
#'
#' @return if the assertion fails an error occurs, otherwise `NULL` is returned
#' invisibly.
#'
#' @examples
#' x = matrix(data = c(0,0.2,0.1,0.5,0.05,0.7), nrow = 2, ncol = 3, byrow = TRUE)
#' colnames(x) = c(12, 34, 42)
#' x
#'
#' assert_hazard_matrix(x)
#'
#' @noRd
assert_hazard_matrix = function(x) {
  assert_matrix(x, any.missing = FALSE, min.rows = 1, min.cols = 1, col.names = "named")
  assert_numeric(as.numeric(colnames(x)), unique = TRUE, sorted = TRUE,
                 any.missing = FALSE, null.ok = FALSE)
  stopifnot(x >= 0)

  invisible(NULL)
}
