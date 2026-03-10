#' Extract time points from a probability matrix or vector
#'
#' Internal helper to consistently obtain the time points associated with
#' a probability matrix or vector.
#'
#' @param x (`numeric()` | `matrix()`)\cr
#'  Probability vector (length = time points) or matrix (rows = observations, columns = time points).
#' @param times (`numeric()` | `NULL`)\cr
#'  Optional vector of time points corresponding to `x`.
#'
#' @return A validated numeric vector of time points.
#'
#' @noRd
#' @keywords internal
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
