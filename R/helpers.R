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
