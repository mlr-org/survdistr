test_that("vec_interp returns x unchanged when eval_times = NULL", {
  x = c(1, 0.8, 0.5)
  times = c(0, 1, 2)
  expect_equal(vec_interp(x, times = times, eval_times = NULL, add_times = FALSE), x)
})

test_that("vec_interp interpolation works", {
  # constant
  x = c(1, 0.8, 0.5)
  times = c(0, 1, 2)
  eval_times = c(0.5, 1.5)
  out = vec_interp(x, times, eval_times, constant = TRUE, type = "surv")
  expect_equal(out, c("0.5" = 1.0, "1.5" = 0.8))

  # linear
  out = vec_interp(x, times, eval_times, constant = FALSE, type = "surv")
  expect_equal(out, c("0.5" = 0.9, "1.5" = 0.65))
})

test_that("vec_interp extrapolation works", {
  # before first time
  x = c(0.9, 0.8)
  times = c(1, 2)
  out = vec_interp(x, times, eval_times = c(0, 0.1), constant = TRUE, type = "surv")
  expect_equal(out, c("0" = 1, "0.1" = 1))
  out = vec_interp(x, times, eval_times = c(0, 0.1), constant = FALSE, type = "surv")
  expect_equal(out, c("0" = 1, "0.1" = 0.99))

  # after last time
  out = vec_interp(x, times, eval_times = c(3, 13), constant = TRUE, type = "surv")
  expect_equal(out, c("3" = 0.8, "13" = 0.8))
  out = vec_interp(x, times, eval_times = c(3, 13), constant = FALSE, type = "surv")
  expect_equal(out, c("3" = 0.7, "13" = 0))
  out = vec_interp(1 - x, times, eval_times = c(3, 13), constant = FALSE, type = "cdf")
  expect_equal(out, c("3" = 0.3, "13" = 1))
})

test_that("vec_interp checks", {
  x = c(1, 1.1) # out of [0,1]
  times = c(0, 1)
  expect_error(vec_interp(x, times, eval_times = 0.5, type = "surv"))

  x = c(1, 0.9, 0.95) # increasing S(t)!
  times = c(0, 1, 2)
  expect_error(vec_interp(x, times, eval_times = 1.5, type = "surv"))

  x = c(0, 0.5, 0.4) # decreasing CDF
  times = c(0, 1, 2)
  expect_error(vec_interp(x, times, eval_times = 1.5, type = "cdf"))
})
