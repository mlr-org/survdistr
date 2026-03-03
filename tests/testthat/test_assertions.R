test_that("assert_prob_matrix works", {
  # no data.frame allowed
  expect_error(assert_prob_matrix(data.frame(1)), "Must be of type")

  # generate a survival matrix
  set.seed(42)
  nrows = 100
  ncols = 500
  x = gen_random_mat(nrows, ncols, "surv")

  # times not specified and no column names
  expect_error(assert_prob_matrix(x), "Time points must be provided")

  # time points specified, but are decreasing
  colnames(x) = ncols:1
  expect_error(assert_prob_matrix(x), "Must be sorted")

  # duplicate time points
  colnames(x) = 1:ncols
  colnames(x)[100] = 99
  expect_error(assert_prob_matrix(x), "Contains duplicated values")

  # time points specified and numeric
  times = 1:ncols
  expect_silent(assert_prob_matrix(x, times))

  # assert_prob works for matrix input
  expect_silent(assert_prob(x, times))

  # time points specified as colnames
  colnames(x) = times
  expect_silent(assert_prob_matrix(x))

  # S(t) >= S(t+1)
  x[100, 489] = 1
  expect_error(assert_prob_matrix(x), "Survival probabilities must be")
  # x is a survival matrix so this fails for CDF/CIF as well
  expect_error(assert_prob_matrix(x, type = "cdf"), "CDF probabilities must be")
  expect_error(assert_prob_matrix(x, type = "cif"), "CIF probabilities must be")

  # checks
  x = matrix(0.9)
  expect_error(assert_prob_matrix(x, times = 0), "survival must equal 1")
  expect_error(assert_prob_matrix(x, times = 0, type = "cdf"), "CDF must equal 0")
  expect_error(assert_prob_matrix(x, times = 0, type = "cif"), "CIF must equal 0")
})

test_that("assert_prob_vec works", {
  # no data.frame allowed
  expect_error(assert_prob_vec(data.frame(1)), "Must be of type")

  # generate a survival vector
  set.seed(42)
  n = 25
  x = sort(runif(n), decreasing = TRUE)

  # times not specified and no names
  expect_error(assert_prob_vec(x), "Time points must be provided")

  # time points specified, but are decreasing
  names(x) = n:1
  expect_error(assert_prob_vec(x), "Must be sorted")

  # duplicate time points
  names(x) = 1:n
  names(x)[25] = 24
  expect_error(assert_prob_vec(x), "Contains duplicated values")

  # time points specified and numeric
  times = 1:n
  expect_silent(assert_prob_vec(x, times))

  # assert_prob works for vector input
  expect_silent(assert_prob(x, times))

  # time points specified as names
  names(x) = times
  expect_silent(assert_prob_vec(x))

  # S(t) >= S(t+1)
  x[24] = x[23] + .01
  expect_error(assert_prob_vec(x), "Survival probabilities must be")

  # checks
  expect_error(assert_prob_vec(-0.21, times = 1), "Probabilities must lie in")
  expect_error(assert_prob_vec(c(0.9, 0.91), times = c(1, 2)), "Survival probabilities must be")
  expect_error(assert_prob_vec(c(0.9, 0.8), times = c(1, 2), type = "cdf"), "CDF probabilities must be")
  expect_error(assert_prob_vec(c(0.9, 0.8), times = c(1, 2), type = "cif"), "CIF probabilities must be")
  expect_error(assert_prob_vec(0.1, times = 0), "survival must equal 1")
  expect_error(assert_prob_vec(0.1, times = 0, type = "cdf"), "CDF must equal 0")
  expect_error(assert_prob_vec(0.1, times = 0, type = "cif"), "CIF must equal 0")
})

test_that("assert_nonneg_matrix works", {
  # Valid matrix with non‑negative entries and proper column names
  x = matrix(c(0, 0.2, 0.1, 0.5, 0.05, 0.7), nrow = 2, byrow = TRUE)
  colnames(x) = c(12, 34, 42)
  expect_silent(assert_nonneg_matrix(x))

  # Matrix with zero values only
  x_zero = matrix(0, nrow = 2, ncol = 3)
  colnames(x_zero) = c(1, 2, 3)
  expect_silent(assert_nonneg_matrix(x_zero))

  # Matrix containing a negative value
  x_neg = matrix(c(0, -0.1, 0.2, 0.3), nrow = 2)
  colnames(x_neg) = c(5, 10)
  expect_error(assert_nonneg_matrix(x_neg), "All entries must be non‑negative")

  # fails for invalid column names
  x = matrix(1:6, nrow = 2)
  # Missing column names
  expect_error(assert_nonneg_matrix(x), "Must have colnames")

  # Column names not increasing
  colnames(x) = c(10, 5, 7)
  expect_error(assert_nonneg_matrix(x), "Must be sorted")
})
