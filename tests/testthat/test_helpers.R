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

test_that("trim_duplicates works", {
  x_vec = c(1, 1, 1, 0.8, 0.8, 0.6)
  names(x_vec) = c(0, 1, 2, 3, 4, 5)
  out = trim_duplicates(x_vec)
  expect_named(out, c("x", "times"))
  expect_equal(out$x, c(1, 0.8, 0.6))
  expect_equal(out$times, c(0, 3, 5))

  x_mat = rbind(
    c(1, 1, 1, 0.8, 0.8, 0.6),
    c(1, 1, 1, 0.7, 0.7, 0.4)
  )
  out = trim_duplicates(x = x_mat, times = 1:6)
  expect_named(out, c("x", "times"))
  expect_equal(out$x, x_mat[, c(1, 4, 6), drop = FALSE])
  expect_equal(out$times, c(1, 4, 6))

  # this cannot be trimmed
  x_mat = rbind(
    c(1, 0.8, 0.8, 0.8),
    c(1, 0.7, 0.6, 0.4)
  )
  out = trim_duplicates(x = x_mat, times = 1:4)
  expect_named(out, c("x", "times"))
  expect_equal(out$x, x_mat)
  expect_equal(out$times, 1:4)

  # higher tolerance should lead to more trimming
  surv = gen_random_mat(nrows = 1, ncols = 1000)[1, ]
  expect_equal(length(trim_duplicates(surv, times = 1:1000, tol = 1e-8)$times), 1000)
  expect_lt(length(trim_duplicates(surv, times = 1:1000, tol = 1e-3)$times), 1000)
  expect_lt(length(trim_duplicates(surv, times = 1:1000, tol = 1e-1)$times), 100)
})
