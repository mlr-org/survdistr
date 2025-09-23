library(checkmate)

gen_mat = function(nrows, ncols, type = "surv") {
  x = matrix(runif(nrows * ncols), nrow = nrows, ncol = ncols)

  # Sort each row in decreasing order => survival matrix
  if (type == "surv") {
    t(apply(x, 1, function(row) sort(row, decreasing = TRUE)))
  } else { # CDF/CIF => increasing values
    t(apply(x, 1, function(row) sort(row, decreasing = FALSE)))
  }
}
