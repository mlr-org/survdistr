test_that("extract_times works", {
  ## derive times from colnames (matrix)
  x_mat = matrix(1:6, nrow = 2)
  colnames(x_mat) = c("0", "1", "2")
  expect_equal(extract_times(x_mat), c(0, 1, 2))

  ## derive times from names (vector)
  x_vec = c(0.9, 0.8, 0.7)
  names(x_vec) = c("0", "1", "2")
  expect_equal(extract_times(x_vec), c(0, 1, 2))

  ## use supplied times (matrix + vector)
  times = c(0, 1, 2)
  expect_equal(extract_times(x_mat, times), times)
  expect_equal(extract_times(x_vec, times), times)

  ## fails if names missing and times NULL
  expect_error(extract_times(c(0.9, 0.8, 0.7)))
  expect_error(extract_times(matrix(1:6, nrow = 2)))

  ## fails on length mismatch
  expect_error(extract_times(x_vec, c(0, 1)))

  ## fails if times not sorted
  expect_error(extract_times(x_vec, c(0, 2, 1)))

  ## fails if times contain duplicates
  expect_error(extract_times(x_vec, c(0, 1, 1)))

  ## fails if times contain negative values
  expect_error(extract_times(x_vec, c(-1, 0, 1)))

  ## fails if x invalid type
  expect_error(extract_times("not a matrix"))
  expect_error(extract_times(list(1, 2, 3)))
})
