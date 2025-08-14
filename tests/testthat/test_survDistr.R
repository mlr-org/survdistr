# TODO Notes:
# new_times values entirely before times.
# new_times values entirely after times.
# Mixed new_times values spanning the range of times.
# Empty times, new_times, or x.
# Check the behavior for both constant and linear interpolation.
set.seed(42)
x = gen_surv_mat(nrows = 2, ncols = 3)
colnames(x) = c(12, 34, 42)
obj = survDistr$new(x)

test_that("constructor", {
  # Valid input
  checkmate::expect_r6(obj, "survDistr")
  expect_equal(obj$data_type, "surv")
  expect_equal(obj$get_data, x)
  expect_equal(obj$times, c(12, 34, 42))
  expect_equal(obj$interp_meth, "const_surv")

  # Invalid inputs
  expect_error(survDistr$new(x = NULL), "Must be of type")
  expect_error(survDistr$new(x = "not a matrix"), "Must be of type")
  expect_error(survDistr$new(x = x[, 0, drop = FALSE]), "Must have at least 1 cols")
  expect_error(survDistr$new(x = matrix()), "Must have colnames")
  expect_error(survDistr$new(x = matrix(dimnames = list(NULL, 1))), "Contains missing values")
})

test_that("constant survival", {
  # Case: `new_times` entirely before `times`
  res = obj$survival(times = 5)
  expect_matrix(res, nrows = 2, ncols = 1, col.names = "named")
  expect_equal(res[, 1], c(1, 1))

  # Case: `new_times` entirely after `times`
  res = obj$survival(times = 50)
  expect_matrix(res, nrows = 2, ncols = 1, col.names = "named")
  expect_equal(res[, 1], obj$get_data[, 3]) # S(t) = S(t_max)

  # Case: Mixed `new_times` (each time point is a bit after the one we have)
  res = obj$survival(times = c(13, 38, 45))
  expect_matrix(res, nrows = 2, ncols = 3, col.names = "named")
  expect_equal(unname(res), obj$data)

  # Case: Empty `new_times`
  res = obj$survival(times = NULL)
  expect_equal(res, obj$get_data) # Returns original matrix
})

test_that("linear survival", {
  obj$interp_meth = "linear_surv"

  # Case: `new_times` entirely before `times`
  res = obj$survival(times = 5)
  expect_matrix(res, nrows = 2, ncols = 1, col.names = "named")
  expect_true(all(res[, 1] < 1))

  # Case: `new_times` entirely after `times`
  res = obj$survival(times = 50)
  expect_matrix(res, nrows = 2, ncols = 1, col.names = "named")
  expect_true(all(res[, 1] < obj$get_data[, 3])) # S(t) < S(t_max)

  # Case: Mixed `new_times`
  res = obj$survival(times = c(13, 38, 45)) # each time is a bit after the one we have
  expect_matrix(res, nrows = 2, ncols = 3, col.names = "named")
  expect_true(all(res < obj$data)) # Linear interpolation within range

  # Case: Empty `new_times`
  res = obj$survival(times = NULL)
  expect_equal(res, obj$get_data) # Returns original matrix
})
