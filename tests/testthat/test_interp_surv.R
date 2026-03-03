# we only test the vector input here since it uses the C matrix version internally
test_that("interp() returns x unchanged when eval_times = NULL", {
  x = c(1, 0.8, 0.5)
  times = c(0, 1, 2)
  expect_equal(interp(x, times = times, eval_times = NULL, add_times = FALSE), x)
})

test_that("interp() checks", {
  x = c(1, 1.1) # out of [0,1]
  times = c(0, 1)
  expect_error(interp(x, times, eval_times = 0.5, type = "surv"))

  x = c(1, 0.9, 0.95) # increasing S(t)!
  times = c(0, 1, 2)
  expect_error(interp(x, times, eval_times = 1.5, type = "surv"))
})

test_that("interp() works", {
  # constant survival
  x = c(1, 0.8, 0.5)
  times = c(0, 1, 2)
  eval_times = c(0.5, 1.5)
  out = interp(x, times, eval_times, method = "const_surv")
  expect_equal(out, c("0.5" = 1.0, "1.5" = 0.8))

  # linear survival
  out = interp(x, times, eval_times, method = "linear_surv")
  expect_equal(out, c("0.5" = 0.9, "1.5" = 0.65))
})

test_that("interp() extrapolation works", {
  # before first time
  x = c(0.9, 0.8)
  times = c(1, 2)
  out = interp(x, times, eval_times = c(0, 0.1), method = "const_surv")
  expect_equal(out, c("0" = 1, "0.1" = 1))
  out = interp(x, times, eval_times = c(0, 0.1), method = "linear_surv")
  expect_equal(out, c("0" = 1, "0.1" = 0.99))

  # after last time
  out = interp(x, times, eval_times = c(3, 13), method = "const_surv")
  expect_equal(out, c("3" = 0.8, "13" = 0.8))
  out = interp(x, times, eval_times = c(3, 13), method = "linear_surv")
  expect_equal(out, c("3" = 0.7, "13" = 0))
})

