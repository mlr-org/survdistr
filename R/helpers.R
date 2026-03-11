#' @title Extract time points from a probability matrix or vector
#'
#' @description
#' Helper function to consistently obtain and validate the time points
#' associated with a probability matrix or vector.
#'
#' @param x (`numeric()` | `matrix()`)\cr
#'  Probability vector (length = time points) or matrix
#'  (rows = observations, columns = time points).
#' @param times (`numeric()` | `NULL`)\cr
#'  Optional vector of time points corresponding to `x`.
#'
#' @return A validated numeric vector of time points.
#' @export
extract_times = function(x, times = NULL) {
  is_mat = is.matrix(x)

  # dimension of time axis
  n_times = if (is_mat) ncol(x) else length(x)
  x_names = if (is_mat) colnames(x) else names(x)

  if (is.null(times)) {
    if (is.null(x_names)) {
      stop("Time points must be provided via 'times' or names/colnames of 'x'.")
    }
    times = assert_numeric(
      as.numeric(x_names),
      lower = 0, unique = TRUE, sorted = TRUE,
      any.missing = FALSE
    )
  } else {
    times = assert_numeric(
      times,
      lower = 0, unique = TRUE, sorted = TRUE,
      any.missing = FALSE, null.ok = FALSE
    )
  }

  if (length(times) != n_times) {
    stop("Length of 'times' must match the time dimension of 'x'.")
  }

  times
}

#' Transform interpolation output
#'
#' Applies the requested output transformation from survival S(t) to other functions
#' (eg. F(t) or H(t)) and optionally attaches time labels to the result.
#'
#' @param res (`numeric` | `matrix`) Interpolated survival values.
#' @param times (`numeric`) Time points corresponding to columns / elements.
#' @param output (`character(1)`) Output type: `"surv"`, `"cdf"`, or `"cumhaz"`.
#' @param add_times (`logical(1)`) Whether to attach time labels.
#' @param eps (`numeric(1)`) Small value used to avoid `-Inf` in `cumhaz`.
#'
#' @return Transformed object with optional time labels.
#'
#' @noRd
#' @keywords internal
transform_result = function(res, times, output, add_times, eps) {
  # transform S(t) output
  if (output == "cdf") {
    res = 1 - res
  } else if (output == "cumhaz") {
    # H(t) = -log(S(t)), floored at eps to avoid -log(0) = Inf
    res = -log(pmax(res, eps))
  }

  # attach time labels
  if (add_times) {
    if (is.matrix(res)) {
      colnames(res) = as.character(times)
    } else {
      names(res) = as.character(times)
    }
  }

  res
}

#' Map interpolation method to internal implementation
#'
#' Maps a user-specified interpolation method to the corresponding internal method.
#' Some methods are aliases for others.
#'
#' @template param_method
#' @return (`character(1)`) Internal method name.
#'
#' @noRd
#' @keywords internal
map_interp_method = function(method) {
  method = assert_choice(method, c("const_surv", "const_dens", "linear_surv", "const_haz", "exp_surv"))

  # keep only the constant aliases
  if (method == "linear_surv") return("const_dens")
  if (method == "exp_surv") return("const_haz")

  method
}

#' @title Remove adjacent duplicate values
#'
#' @description
#' Removes adjacent duplicate values over the time axis, possibly from a
#' probability vector or matrix (e.g. survival curves).
#' Equality is determined with a numeric tolerance.
#'
#' For matrices, duplicate detection is done column-wise across all rows.
#' Only the earliest time point in each run of (near-)equal values is kept.
#'
#' @param x (`numeric()` | `matrix()`)
#'  Vector (length = time points) or matrix (rows = observations, columns = time points).
#' @param times (`numeric()` | `NULL`)
#'  Optional time points corresponding to `x`.
#'  If `NULL`, extracted from names/colnames.
#' @param tol (`numeric(1)`)
#'  Absolute tolerance used to detect equality between adjacent time points.
#'
#' @return A named list with:
#'  * `x`: numeric vector or matrix with duplicate adjacent time points removed.
#'  * `times`: numeric vector of retained time points.
#' @examples
#' # remove adjacent duplicates from a survival probability vector
#' surv = c(1, 1, 0.8, 0.8, 0.5, 0.5, 0.2)
#' trim_duplicates(surv, times = 1:7)
#'
#' @export
trim_duplicates = function(x, times = NULL, tol = 1e-10) {
  is_mat = is.matrix(x)
  if (is_mat) {
    assert_matrix(x, mode = "numeric", any.missing = FALSE, min.rows = 1, min.cols = 1)
  } else {
    x = assert_numeric(x, any.missing = FALSE, min.len = 1)
  }

  tol = assert_numeric(tol, lower = 0, finite = TRUE, len = 1)
  times = extract_times(x, times)

  # remove times
  if (is_mat) {
    colnames(x) = NULL
  } else {
    names(x) = NULL
  }

  n_times = length(times)
  if (n_times == 1L) {
    return(list(x = x, times = times))
  }

  keep = rep(FALSE, n_times)
  keep[1] = TRUE
  ref_idx = 1L

  for (j in 2:n_times) {
    is_dup = if (is_mat) {
      all(abs(x[, j] - x[, ref_idx]) <= tol)
    } else {
      abs(x[j] - x[ref_idx]) <= tol
    }

    if (!is_dup) {
      keep[j] = TRUE
      ref_idx = j
    }
  }

  if (is_mat) {
    x = x[, keep, drop = FALSE]
  } else {
    x = x[keep]
  }
  times = times[keep]

  list(x = x, times = times)
}
