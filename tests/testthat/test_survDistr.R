# TODO Notes:
# new_times values entirely before times.
# new_times values entirely after times.
# Mixed new_times values spanning the range of times.
# Empty times, new_times, or x.
# Check the behavior for both constant and linear interpolation.
set.seed(42)
x = gen_mat(nrows = 2, ncols = 3, type = "surv")
times = c(12, 34, 42)
obj = survDistr$new(x, times)

test_that("constructor", {
  # Valid input
  checkmate::expect_r6(obj, "survDistr")
  expect_equal(obj$data_type, "surv")
  expect_equal(obj$data(add_times = FALSE), x)
  expect_equal(obj$times, c(12, 34, 42))
  expect_equal(obj$interp_meth, "const_surv")

  # Invalid inputs
  expect_error(survDistr$new(x = NULL), "Must be of type")
  expect_error(survDistr$new(x = "not a matrix"), "Must be of type")
  expect_error(survDistr$new(x = x[, 0, drop = FALSE]), "Must have at least 1 cols")
  expect_error(survDistr$new(x = matrix()), "Contains missing values")
  expect_error(survDistr$new(x = matrix(dimnames = list(NULL, 1))), "Contains missing values")
})

test_that("constant survival", {
  # Extrapolation case 1: requested time before earliest time
  res = obj$survival(times = 5)
  expect_matrix(res, nrows = 2, ncols = 1, col.names = "named")
  expect_equal(res[, 1], c(1, 1)) # S(t < t_min) = 1

  # Extrapolation case 2: requested time after latest time
  res = obj$survival(times = 50)
  expect_matrix(res, nrows = 2, ncols = 1, col.names = "named")
  expect_equal(res[, 1], obj$data()[, 3]) # S(t > t_max) = S(t_max)

  # Interpolation case: requested time points are a bit after the original ones we have
  res = obj$survival(times = c(13, 38, 45), add_times = FALSE)
  expect_matrix(res, nrows = 2, ncols = 3, col.names = "unnamed")
  expect_equal(res, obj$data(add_times = FALSE)) # equal due to left-continuous interpolation

  # Case: Empty `new_times`
  res = obj$survival(times = NULL)
  expect_equal(res, obj$data()) # Returns original matrix
})

test_that("linear survival", {
  # can 'hack' the interpolation method
  obj$interp_meth = "linear_surv"

  # Extrapolation case 1: requested time before earliest time
  res = obj$survival(times = 5)
  expect_matrix(res, nrows = 2, ncols = 1, col.names = "named")
  expect_true(all(res[, 1] < 1)) # S(t > t_max) = 1

  # Extrapolation case 2: requested time after latest time
  res = obj$survival(times = 50)
  expect_matrix(res, nrows = 2, ncols = 1, col.names = "named")
  expect_true(all(res[, 1] < obj$data()[, 3])) # S(t > t_max) < S(t_max)

  # Interpolation case: requested time points are a bit after the original ones we have
  res = obj$survival(times = c(13, 38, 45))
  expect_matrix(res, nrows = 2, ncols = 3, col.names = "named")
  expect_true(all(res < obj$data()))

  # Case: Empty `new_times`
  res = obj$survival(times = NULL)
  expect_equal(res, obj$data()) # Returns original matrix
})
