test_that("matrix: derive times from colnames", {
  x = matrix(1:6, nrow = 2)
  colnames(x) = c("0", "1", "2")

  times = extract_times(x)

  expect_equal(times, c(0, 1, 2))
})

test_that("vector: derive times from names", {
  x = c(0.9, 0.8, 0.7)
  names(x) = c("0", "1", "2")

  times = extract_times(x)

  expect_equal(times, c(0, 1, 2))
})

test_that("uses supplied times (matrix)", {
  x = matrix(1:6, nrow = 2)
  times = c(0, 1, 2)

  expect_equal(extract_times(x, times), times)
})

test_that("uses supplied times (vector)", {
  x = c(0.9, 0.8, 0.7)
  times = c(0, 1, 2)

  expect_equal(extract_times(x, times), times)
})

test_that("fails if names missing and times NULL", {
  x_vec = c(0.9, 0.8, 0.7)
  x_mat = matrix(1:6, nrow = 2)

  expect_error(extract_times(x_vec))
  expect_error(extract_times(x_mat))
})

test_that("fails on length mismatch", {
  x_vec = c(0.9, 0.8, 0.7)
  names(x_vec) = c("0", "1", "2")

  expect_error(extract_times(x_vec, c(0, 1)))
})

test_that("fails if times not sorted", {
  x = c(0.9, 0.8, 0.7)
  expect_error(extract_times(x, c(0, 2, 1)))
})

test_that("fails if times contain duplicates", {
  x = c(0.9, 0.8, 0.7)
  expect_error(extract_times(x, c(0, 1, 1)))
})

test_that("fails if times contain negative values", {
  x = c(0.9, 0.8, 0.7)
  expect_error(extract_times(x, c(-1, 0, 1)))
})

test_that("fails if x is not numeric or matrix", {
  expect_error(extract_times("not a matrix"))
  expect_error(extract_times(list(1, 2, 3)))
})
