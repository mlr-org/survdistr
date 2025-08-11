# If `x` is a matrix with monotonically non-decreasing values in each row,
# e.g. a CDF matrix, then this function calculates the PDF.
# i.e. this is a F(t) => f(t) transformation for the discrete time setting.
# f(t_k) = F(t_k) - F(t_{k-1})
.rowwise_diffs = function(x) {
  d = x[, -1, drop = FALSE] - x[, -ncol(x), drop = FALSE]
  # Add the first column
  cbind(x[, 1, drop = FALSE], d)
}
