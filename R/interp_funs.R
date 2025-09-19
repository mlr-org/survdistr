#' Interpolate Survival / CDF / CIF Matrices
#'
#' Wrapper around the internal C++ interpolation function \code{c_mat_interp}.
#' Performs input validation before calling the underlying C++ code.
#' Can be used for survival, cumulative distribution (CDF), or cumulative
#' incidence (CIF) matrices.
#'
#' @param x (`matrix()`)\cr Survival/CDF/CIF matrix with rows as observations
#'  and columns as time points.
#' @param times (`numeric()`|`NULL`)\cr
#'  Original time points corresponding to columns of `x`.
#' @param eval_times (`numeric()`|`NULL`)\cr
#'  New time points at which to interpolate.
#'  Values do not need to be sorted or unique, just non-negative.
#'  If `NULL`, `x` is returned unchanged.
#' @param constant (`logical(1)`)\cr
#'  If `TRUE`, use piecewise-constant (left-continuous) interpolation.
#'  If `FALSE`, use piecewise-linear interpolation.
#' @param type (`character(1)`)\cr
#'  One of `"surv"`, `"cdf"`, or `"cif"`, indicating the data type of the input matrix.
#' @param add_times (`logical(1)`)\cr
#'  If `TRUE`, sets the column names of the returned matrix to the requested `eval_times`.
#' @param check (`logical(1)`)\cr
#'  If `TRUE`, run input matrix validation via [assert_surv_matrix()];
#'  set `FALSE` to skip checks (NOT recommended for external use).
#'
#' @return A numeric matrix with the same number of rows as `x` and number of
#'   columns equal to `length(eval_times)`.
#'
#' @examples
#' x = matrix(c(1,0.8,0.6,
#'              1,0.7,0.4), nrow = 2, byrow = TRUE)
#' times = c(0, 10, 20)
#' eval_times = c(5, 15, 25, 15)  # duplicates & unordered
#' mat_interp(x, times, eval_times, constant = TRUE, type = "surv")
#' @export
mat_interp = function(x, times = NULL, eval_times = NULL, constant = TRUE, type = "surv",
                      add_times = TRUE, check = TRUE) {
  # general checks
  eval_times = assert_numeric(eval_times, lower = 0, any.missing = FALSE,
                              null.ok = TRUE, min.len = 1)
  if (is.null(eval_times)) {
    return(x)
  }
  assert_flag(constant)
  type = assert_choice(type, c("surv", "cdf", "cif"))
  assert_flag(add_times)
  assert_flag(check)

  # Optional matrix check
  if (check) {
    if (type == "surv") {
      times = assert_surv_matrix(x, times)
    } else {
      # Transform cdf/cif to survival for checks
      times = assert_surv_matrix(1 - x, times)
    }
  }

  # unique + sorted eval_times for C++
  eval_times_unique = sort(unique(eval_times))

  # call C++ interpolation once
  mat = c_mat_interp(x, times, eval_times_unique, constant, type)

  # map back to requested order (with duplicates) if necessary
  if (!identical(eval_times, eval_times_unique)) {
    idx = match(eval_times, eval_times_unique)
    mat = mat[, idx, drop = FALSE]
  }

  if (add_times) {
    colnames(mat) = eval_times
  }

  mat
}
