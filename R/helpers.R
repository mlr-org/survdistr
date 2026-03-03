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
