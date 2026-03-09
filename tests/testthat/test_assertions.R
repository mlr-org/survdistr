test_that("assert_prob works with input matrix", {
  # no data.frame allowed
  expect_error(assert_prob(data.frame(1)), "Must be of type")

  # generate a survival matrix
  set.seed(42)
  nrows = 100
  ncols = 500
  x = gen_random_mat(nrows, ncols, "surv")

  # times not specified and no column names
  expect_error(assert_prob(x), "Time points must be provided")

  # time points specified, but are decreasing
  colnames(x) = ncols:1
  expect_error(assert_prob(x), "Must be sorted")

  # duplicate time points
  colnames(x) = 1:ncols
  colnames(x)[100] = 99
  expect_error(assert_prob(x), "Contains duplicated values")

  # time points specified and numeric
  times = 1:ncols
  expect_silent(assert_prob(x, times))

  # assert_prob works for matrix input
  expect_silent(assert_prob(x, times))

  # time points specified as colnames
  colnames(x) = times
  expect_silent(assert_prob(x))

  # S(t) >= S(t+1)
  x[100, 489] = 1
  expect_error(assert_prob(x), "Survival probabilities must be")
  # x is a survival matrix so this fails for CDF/CIF as well
  expect_error(assert_prob(x, type = "cdf"), "CDF probabilities must be")
  expect_error(assert_prob(x, type = "cif"), "CIF probabilities must be")

  # checks
  x = matrix(0.9)
  expect_error(assert_prob(x, times = 0), "survival must equal 1")
  expect_error(assert_prob(x, times = 0, type = "cdf"), "CDF must equal 0")
  expect_error(assert_prob(x, times = 0, type = "cif"), "CIF must equal 0")
  expect_error(assert_prob(x, times = 0, type = "haz"), "times must be strictly positive.")
  expect_error(assert_prob(x, times = 0, type = "dens"), "times must be strictly positive.")
})

test_that("assert_prob works with input vector", {
  # no data.frame allowed
  expect_error(assert_prob(data.frame(1)), "Must be of type")

  # generate a survival vector
  set.seed(42)
  n = 25
  x = sort(runif(n), decreasing = TRUE)

  # times not specified and no names
  expect_error(assert_prob(x), "Time points must be provided")

  # time points specified, but are decreasing
  names(x) = n:1
  expect_error(assert_prob(x), "Must be sorted")

  # duplicate time points
  names(x) = 1:n
  names(x)[25] = 24
  expect_error(assert_prob(x), "Contains duplicated values")

  # time points specified and numeric
  times = 1:n
  expect_silent(assert_prob(x, times))

  # assert_prob works for vector input
  expect_silent(assert_prob(x, times))

  # time points specified as names
  names(x) = times
  expect_silent(assert_prob(x))

  # S(t) >= S(t+1)
  x[24] = x[23] + .01
  expect_error(assert_prob(x), "Survival probabilities must be")

  # checks
  expect_error(assert_prob(-0.21, times = 1), "Survival probabilities must be")
  expect_error(assert_prob(c(0.9, 0.91), times = c(1, 2)), "Survival probabilities must be")
  expect_error(assert_prob(c(0.9, 0.8), times = c(1, 2), type = "cdf"), "CDF probabilities must be")
  expect_error(assert_prob(c(0.9, 0.8), times = c(1, 2), type = "cif"), "CIF probabilities must be")
  expect_error(assert_prob(0.1, times = 0), "survival must equal 1")
  expect_error(assert_prob(0.1, times = 0, type = "cdf"), "CDF must equal 0")
  expect_error(assert_prob(0.1, times = 0, type = "cif"), "CIF must equal 0")
  expect_error(assert_prob(1.2, times = 1, type = "haz"),
               "Hazard probabilities must be in [0,1].", fixed = TRUE)
  expect_error(assert_prob(-0.1, times = 1, type = "haz"),
               "Hazard probabilities must be in [0,1].", fixed = TRUE)
  expect_error(assert_prob(1.2, times = 1, type = "dens"),
               "Density probabilities must be in [0,1].", fixed = TRUE)
  expect_error(assert_prob(-0.1, times = 1, type = "dens"),
               "Density probabilities must be in [0,1].", fixed = TRUE)
})
