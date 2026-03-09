library(checkmate)

gen_random_mat = function(nrows, ncols, type = "surv") {
  x = matrix(runif(nrows * ncols), nrow = nrows, ncol = ncols)

  if (type == "surv") {
    # sort each row in decreasing order => S(t)
    t(apply(x, 1, function(row) sort(row, decreasing = TRUE)))
  } else if (type == "cdf" || type == "cif") {
    # CDF/CIF => increasing values
    t(apply(x, 1, function(row) sort(row, decreasing = FALSE)))
  } else if (type = "prob" || type == "haz" || type == "dens") {
    # no monotonicity requirement, just probabilities
    x
  } else {
    stop("Invalid type")
  }
}
