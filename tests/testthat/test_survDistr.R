set.seed(42)
x = gen_random_mat(nrows = 3, ncols = 3, type = "surv")
times = c(12, 34, 42)
obj = survDistr$new(x, times)

test_that("new() works", {
  # Valid input
  checkmate::expect_r6(obj, "survDistr")
  expect_equal(obj$data(add_times = FALSE), x)
  expect_equal(obj$times, c(12, 34, 42))
  expect_equal(obj$interp_meth, "const_surv")
  expect_silent(survDistr$new(x, times = times, check = FALSE)) # skip checks, should still work

  # Invalid inputs
  expect_error(survDistr$new(x = NULL), "Must be of type")
  expect_error(survDistr$new(x = "not a matrix"), "Must be of type")
  expect_error(survDistr$new(x = x[, 0, drop = FALSE]), "Must have at least 1 cols")
  expect_error(survDistr$new(x = matrix()), "Contains missing values")
  expect_error(survDistr$new(x = matrix(dimnames = list(NULL, 1))), "Contains missing values")
})

test_that("print() works", {
  expect_output(print(obj), "survival matrix")
})

test_that("filter() works", {
  obj2 = obj$clone(deep = TRUE)

  # can't filter out of bounds (2 observations only)
  expect_error(obj2$filter(rows = c(0, 2)), ">= 1")
  expect_error(obj2$filter(rows = c(1, 4)), "<= 3")
  expect_error(obj2$filter(rows = c(FALSE, TRUE)), "Must have length 3")
  expect_error(obj2$filter(rows = c(1, 1)), "duplicated values")
  expect_error(obj2$filter(rows = c(2, 1)), "be sorted")

  # no filtering => same data
  expect_equal(obj2$filter()$data(), obj2$data())

  # filter to 2 observations
  expect_invisible(obj2$filter(rows = c(1, 3)))
  expect_equal(obj2$data(add_times = FALSE), x[c(1, 3), , drop = FALSE])
  # filter to 1 observation
  obj2$filter(rows = c(FALSE, TRUE))
  expect_equal(obj2$data(add_times = FALSE), x[3, , drop = FALSE])
  # remove last observation
  obj2$filter(rows = FALSE)
  expect_equal(dim(obj2$data(add_times = FALSE)), c(0, 3))
  # obj remains unchanged
  expect_equal(obj$data(add_times = FALSE), x)
})

test_that("subsetting using 'rows' works", {
  obj3 = obj$clone(deep = TRUE)

  expect_equal(obj3$data(rows = 1, add_times = FALSE), x[1, , drop = FALSE])
  expect_equal(obj3$survival(rows = c(1, 3), add_times = FALSE), x[c(1, 3), , drop = FALSE])
  expect_equal(obj3$cdf(rows = c(1, 3), add_times = FALSE), 1 - x[c(1, 3), , drop = FALSE])
  # obj3 remains unchanged
  expect_equal(obj3$data(add_times = FALSE), x)
})

test_that("survival() works", {
  obj2 = obj$clone(deep = TRUE)

  # constant survival interpolation (default)
  t = c(0, 7, 12, 22, 34, 40, 42, 50)
  res = obj2$survival(times = t)
  expect_matrix(res, nrows = 3, ncols = length(t), col.names = "named")
  res2 = obj2$survival()
  expect_equal(res2, obj2$data()) # Returns original matrix

  # linear survival interpolation
  obj2$interp_meth = "linear_surv"
  res3 = obj2$survival(times = t)
  time_cols = as.character("0", "12", "34", "42")
  expect_equal(res3[, time_cols], res[, time_cols]) # no interpolation at original time points
})

test_that("cdf() works", {
  t = c(0, 7, 12, 22, 34, 40, 42, 50)
  res = obj$cdf(times = t)
  expect_matrix(res, nrows = 3, ncols = length(t), col.names = "named")

  # Returns original matrix transformed to CDF
  res2 = obj$cdf()
  expect_equal(res2, 1 - obj$data())
})

test_that("cumhazard() works", {
  obj2 = obj$clone(deep = TRUE)
  obj2$interp_meth = "linear_surv" # H(t) grows

  t = c(0, 7, 12, 22, 34, 40, 42, 50)
  res = obj2$cumhazard(times = t)
  expect_matrix(res, nrows = 3, ncols = length(t), col.names = "named")
  expect_equal(res[,1], c(0, 0, 0)) # H(0) = 0
  expect_all_true(res[1, ] >= 0) # H(t) increases
  expect_all_true(res[2, ] >= 0)
  expect_all_true(res[3, ] >= 0)

  # large times so that eps kicks in
  res2 = obj2$cumhazard(times = 1000)
  res3 = obj2$cumhazard(times = 1000, eps = 1e-10) # lower eps => lower S(t) => higher H(t)
  expect_all_true(res2[,1] < res3[,1])
})
