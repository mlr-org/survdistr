test_that("as_survDistr() works", {
  mat = matrix(c(1, 0.8, 0.5, 1, 0.9, 0.7), nrow = 2, byrow = TRUE)
  obj = as_survDistr(mat, times = c(1, 2, 3))

  checkmate::expect_r6(obj, "survDistr")
  expect_equal(obj$times, c(1, 2, 3))
  expect_equal(obj$data(add_times = FALSE), mat)

  # keeps existing survDistr objects unchanged
  mat = matrix(c(1, 0.8, 0.5), nrow = 1)
  x = survDistr$new(mat, times = c(1, 2, 3))
  expect_identical(as_survDistr(x), x)

  # invalid input
  expect_error(as_survDistr(1:3), "No as_survDistr() method for objects of class", fixed = TRUE)
})
