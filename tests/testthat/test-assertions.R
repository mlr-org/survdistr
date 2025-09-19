test_that("assert_surv_matrix works", {
  # no data.frame allowed
  expect_error(assert_surv_matrix(data.frame(1)), "Must be of type")

  # make a survival matrix
  set.seed(42)
  nrows = 100
  ncols = 500
  x = gen_surv_mat(nrows, ncols)

  # times not specified and no column names
  expect_error(assert_surv_matrix(x), "Column names are required")

  # time points specified, but are decreasing
  colnames(x) = ncols:1
  expect_error(assert_surv_matrix(x), "Must be sorted")

  # duplicate time points
  colnames(x) = 1:ncols
  colnames(x)[100] = 99
  expect_error(assert_surv_matrix(x), "Contains duplicated values")

  # time points specified and numeric
  times = 1:ncols
  expect_silent(assert_surv_matrix(x, times))

  # time points specified as colnames
  colnames(x) = times
  expect_silent(assert_surv_matrix(x))

  # S(t) >= S(t+1)
  x[100, 489] = 1
  expect_error(assert_surv_matrix(x), "Survival probabilities must be")
})
