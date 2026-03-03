# we only test the vector input here since it uses the C matrix version internally
test_that("interp() returns x unchanged when eval_times = NULL", {
  x = c(1, 0.8, 0.5)
  times = c(0, 1, 2)
  expect_equal(interp(x, times = times, eval_times = NULL, add_times = FALSE), x)
  expect_named(interp(x, times = times, eval_times = NULL, add_times = TRUE))
})

test_that("interp() checks", {
  x = c(1, 1.1) # out of [0,1]
  times = c(0, 1)
  expect_error(interp(x, times, eval_times = 0.5))

  x = c(1, 0.9, 0.95) # increasing S(t)!
  times = c(0, 1, 2)
  expect_error(interp(x, times, eval_times = 1.5))
})

test_that("interp() S(t) works", {
  # S(t) all different, 3 time points, t1 == 0
  x = c(1, 0.8, 0.6)
  times = c(0, 1, 2)
  eval_times = c(0, 0.5, 1, 1.5, 2, 3)
  out = interp(x, times, eval_times, method = "const_surv", add_times = FALSE)
  expect_equal(out, c(1.0, 1.0, 0.8, 0.8, 0.6, 0.6))

  eval_times = c(0, 0.5, 1, 1.5, 2, 3, 4, 5, 6)
  out = interp(x, times, eval_times, method = "linear_surv", add_times = FALSE)
  expect_equal(out, c(1.0, 0.9, 0.8, 0.7, 0.6, 0.4, 0.2, 0.0, 0.0))

  # output CDF and H(t) are transformed correctly
  out_cdf = interp(x, times, eval_times, method = "linear_surv", output = "cdf", add_times = FALSE)
  expect_equal(out_cdf, 1 - out)
  out_cumhaz = interp(x, times, eval_times, method = "linear_surv", output = "cumhaz", add_times = FALSE)
  expect_equal(out_cumhaz, -log(pmax(out, 1e-6)))

  # S(t) all different, 3 time points, t1 > 0
  x = c(0.9, 0.7, 0.5) # slope = 0.2
  times = c(1, 2, 3)
  eval_times = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6)
  out = interp(x, times, eval_times, method = "const_surv", add_times = FALSE)
  expect_equal(out, c(1.0, 1.0, 0.9, 0.9, 0.7, 0.7, 0.5, 0.5, 0.5, 0.5))
  out = interp(x, times, eval_times, method = "linear_surv", add_times = FALSE)
  expect_equal(out, c(1.0, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.3, 0.1, 0.0))

  # S(t) all different, 2 time points, t1 > 0
  x = c(0.8, 0.3) # slope = 0.25
  times = c(1, 2)
  eval_times = c(0, 0.5, 1, 1.5, 2, 2.5, 3)
  out = interp(x, times, eval_times, method = "const_surv", add_times = FALSE)
  expect_equal(out, c(1.0, 1.0, 0.8, 0.8, 0.3, 0.3, 0.3))
  out = interp(x, times, eval_times, method = "linear_surv", add_times = FALSE)
  expect_equal(out, c(1.0, 0.9, 0.8, 0.55, 0.3, 0.05, 0.0))

  # 1 time point only, at t = 0
  out = interp(x = 1, times = 0, eval_times = c(0, 1), method = "const_surv", add_times = FALSE)
  expect_equal(out, c(1, 1))
  out = interp(x = 1, times = 0, eval_times = c(0, 1), method = "linear_surv", add_times = FALSE)
  expect_equal(out, c(1, 1))

  # 1 time point only, at t != 0
  out = interp(x = 0.8, times = 1, eval_times = c(0, 0.5, 1, 1.5, 2),
               method = "const_surv", add_times = FALSE)
  expect_equal(out, c(1, 1, 0.8, 0.8, 0.8))
  out = interp(x = 0.8, times = 1, eval_times = c(0, 0.5, 1, 1.5, 2, 3, 4, 5),
               method = "linear_surv", add_times = FALSE)
  expect_equal(out, c(1, 0.9, 0.8, 0.7, 0.6, 0.4, 0.2, 0.0))
})
