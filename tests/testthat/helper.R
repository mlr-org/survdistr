gen_surv_mat = function(nrows, ncols) {
  x = matrix(runif(nrows * ncols), nrow = nrows, ncol = ncols)

  # Sort each row in decreasing order to make it a survival matrix
  t(apply(x, 1, function(row) sort(row, decreasing = TRUE)))
}
