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

#' Row-wise first differences of a matrix
#'
#' If `x` is a matrix with monotonically non-decreasing values in each row
#' (e.g., a CDF matrix), this function calculates the corresponding PDF
#' by computing row-wise first differences.
#'
#' This corresponds to the transformation \( F(t) \;\Rightarrow\; f(t) \)
#' in the discrete-time setting:
#' \deqn{f(t_k) = F(t_k) - F(t_{k-1})}{f(t_k) = F(t_k) - F(t_{k-1})}
#'
#' @param x A numeric matrix with non-decreasing values in each row.
#'
#' @return A numeric matrix of the same dimensions as `x`, containing the
#'   row-wise differences (PDF values).
#'
#' @examples
#' # Example: simple CDF matrix
#' cdf_mat = matrix(c(0.2, 0.5, 0.8, 0.1, 0.4, 0.9), nrow = 2, byrow = TRUE)
#' cdf_mat
#'
#' # Compute the corresponding PDF matrix
#' rowwise_diffs(cdf_mat)
#'
#' @noRd
#' @keywords internal
rowwise_diffs = function(x) {
  d = x[, -1, drop = FALSE] - x[, -ncol(x), drop = FALSE]

  # Add the first column
  cbind(x[, 1, drop = FALSE], d)
}
