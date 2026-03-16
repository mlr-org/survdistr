# we only test the vector input here since it uses the C matrix version internally
test_that("interp() density works", {
  # S(t) all different, 3 time points, t1 == 0
  x = c(1, 0.8, 0.5)
  times = c(0, 1, 2)
  eval_times = c(0, 0.5, 1, 1.5, 2, 3)

  out = interp(x, times, eval_times, method = "const_surv", output = "density", add_times = FALSE)
  expect_equal(out, c(0.0, 0.0, 0.2, 0.0, 0.3, 0.0))

  out = interp(x, times, eval_times, method = "const_dens", output = "density", add_times = FALSE)
  expect_equal(out, c(0.2, 0.2, 0.2, 0.3, 0.3, 0.3))
  # S(t) != 0 at the last time point, so we get a non-zero density there
  # f(t=0) follows the first interval value

  out = interp(x, times, eval_times, method = "const_haz", output = "density", add_times = FALSE)
  lambda_1 = -log(0.8/1) * 1/(1-0)
  lambda_2 = -log(0.5/0.8) * 1/(2-1)
  expect_equal(out,
    c(lambda_1, # t = 0
      lambda_1 * 1 * (0.8/1)^(0.5-0)/(1-0), # t = 0.5
      lambda_1 * 0.8, # t = 1
      lambda_2 * 0.8 * (0.5/0.8)^(1.5-1)/(2-1), # t = 1.5
      lambda_2 * 0.5, # t = 2
      lambda_2 * 0.5 * (0.5/0.8)^(3-2)/(2-1) # t = 3
    ), tolerance = 1e-6)

  # S(t) all different, 3 time points, t1 > 0
  x = c(0.9, 0.7, 0.4)
  times = c(1, 2, 3)
  eval_times = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5)

  out = interp(x, times, eval_times, method = "const_surv", output = "density", add_times = FALSE)
  expect_equal(out, c(0.0, 0.0, 0.1, 0.0, 0.2, 0.0, 0.3, 0.0, 0.0))

  out = interp(x, times, eval_times, method = "const_dens", output = "density", add_times = FALSE)
  expect_equal(out, c(0.1, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.3, 0.0)) # S(4) > 0, S(5) = 0

  out = interp(x, times, eval_times, method = "const_haz", output = "density", add_times = FALSE)
  lambda_1 = -log(0.9/1) * 1/(1-0)
  lambda_2 = -log(0.7/0.9) * 1/(2-1)
  lambda_3 = -log(0.4/0.7) * 1/(3-2)
  expect_equal(out,
    c(lambda_1, # t = 0
      lambda_1 * 1 * (0.9/1)^(0.5-0)/(1-0), # t = 0.5
      lambda_1 * 0.9, # t = 1
      lambda_2 * 0.9 * (0.7/0.9)^(1.5-1)/(2-1), # t = 1.5
      lambda_2 * 0.7, # t = 2
      lambda_3 * 0.7 * (0.4/0.7)^(2.5-2)/(3-2), # t = 2.5
      lambda_3 * 0.4, # t = 3
      lambda_3 * 0.4 * (0.4/0.7)^(4-3)/(3-2), # t = 4
      lambda_3 * 0.4 * (0.4/0.7)^(5-3)/(3-2) # t = 5
    ), tolerance = 1e-6)

  # S(t) all different, 2 time points, t1 > 0
  x = c(0.8, 0.3)
  times = c(1, 2)
  eval_times = c(0, 0.5, 1, 1.5, 2, 2.5, 3)

  out = interp(x, times, eval_times, method = "const_surv", output = "density", add_times = FALSE)
  expect_equal(out, c(0.0, 0.0, 0.2, 0.0, 0.5, 0.0, 0.0))

  out = interp(x, times, eval_times, method = "const_dens", output = "density", add_times = FALSE)
  expect_equal(out, c(0.2, 0.2, 0.2, 0.5, 0.5, 0.5, 0.0))

  out = interp(x, times, eval_times, method = "const_haz", output = "density", add_times = FALSE)
  lambda_1 = -log(0.8/1) * 1/(1-0)
  lambda_2 = -log(0.3/0.8) * 1/(2-1)
  expect_equal(out,
    c(lambda_1, # t = 0
      lambda_1 * 1 * (0.8/1)^(0.5-0)/(1-0), # t = 0.5
      lambda_1 * 0.8, # t = 1
      lambda_2 * 0.8 * (0.3/0.8)^(1.5-1)/(2-1), # t = 1.5
      lambda_2 * 0.3, # t = 2
      lambda_2 * 0.3 * (0.3/0.8)^(2.5-2)/(2-1), # t = 2.5
      lambda_2 * 0.3 * (0.3/0.8)^(3-2)/(2-1) # t = 3
  ), tolerance = 1e-6)

  # 1 time point only, at t = 0
  out = interp(x = 1, times = 0, eval_times = c(0, 1), method = "const_surv",
               output = "density", add_times = FALSE)
  expect_equal(out, c(0, 0))

  out = interp(x = 1, times = 0, eval_times = c(0, 1), method = "const_dens",
               output = "density", add_times = FALSE)
  expect_equal(out, c(0, 0))

  out = interp(x = 1, times = 0, eval_times = c(0, 1), method = "const_haz",
               output = "density", add_times = FALSE)
  expect_equal(out, c(0, 0))

  # 1 time point > 0, with S > 0
  out = interp(x = 0.4, times = 1, eval_times = c(0, 0.5, 1, 1.5, 2),
               method = "const_surv", output = "density", add_times = FALSE)
  expect_equal(out, c(0, 0, 0.6, 0, 0))

  out = interp(x = 0.4, times = 1, eval_times = c(0, 0.5, 1, 1.5, 2),
               method = "const_dens", output = "density", add_times = FALSE)
  expect_equal(out, c(0.6, 0.6, 0.6, 0.6, 0.0))

  out = interp(x = 0.4, times = 1, eval_times = c(0, 0.5, 1, 1.5, 2, 10, 20),
               method = "const_haz", output = "density", add_times = FALSE)
  lambda_1 = -log(0.4) # delta_1 = 1
  expect_equal(out,
    c(lambda_1, # t = 0
      lambda_1 * 0.4 ^ (0.5 - 0), # t = 0.5
      lambda_1 * 0.4, # t = 1
      lambda_1 * 0.4 * 0.4^(1.5 - 1), # t = 1.5
      lambda_1 * 0.4 * 0.4^(2 - 1), # t = 2
      lambda_1 * 0.4 * 0.4^(10 - 1), # t = 10
      lambda_1 * 0.4 * 0.4^(20 - 1) # t = 20
    ), tolerance = 1e-6)

  # 1 time point > 0, with S = 0 (super edge case)
  out = interp(x = 0, times = 1, eval_times = c(0, 0.5, 1, 1.5),
               method = "const_surv", output = "density", add_times = FALSE)
  expect_equal(out, c(0, 0, 1.0, 0.0))

  out = interp(x = 0, times = 1, eval_times = c(0, 0.5, 1, 1.5),
               method = "const_dens", output = "density", add_times = FALSE)
  expect_equal(out, c(1.0, 1.0, 1.0, 0.0))

  out = interp(x = 0, times = 1, eval_times = c(0, 0.5, 1, 1.5),
               method = "const_haz", output = "density", add_times = FALSE)
  # S(t) = 0 is replaced by eps for the constant hazard interpolation
  eps = 1e-12
  lambda_1 = -log(eps) # delta_1 = 1
  expect_equal(out, c(lambda_1, lambda_1 * 1 * (eps/1)^0.5, lambda_1 * 0.0, 0.0))

  # S(t) with tied values
  x = c(1, 0.8, 0.8, 0.6)
  times = 1:4
  eval_times = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6)

  out = interp(x, times, eval_times, method = "const_surv", output = "density", add_times = FALSE)
  expect_equal(out, c(0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0))

  out = interp(x, times, eval_times, method = "const_dens", output = "density", add_times = FALSE)
  expect_equal(out, c(0.0, 0.0, 0.0, 0.2, 0.2, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2))

  out = interp(x, times, eval_times, method = "const_haz", output = "density", add_times = FALSE)
  # all deltas are 1
  lambda_1 = -log(1/1)
  lambda_2 = -log(0.8/1)
  lambda_3 = -log(0.8/0.8)
  lambda_4 = -log(0.6/0.8)
  expect_equal(out,
    c(lambda_1, # t = 0
      lambda_1 * 1 * (1/1)^(0.5-0), # t = 0.5
      lambda_1 * 1, # t = 1
      lambda_2 * 1 * (0.8/1)^(1.5-1), # t = 1.5
      lambda_2 * 0.8, # t = 2
      lambda_3 * 0.8 * (0.8/0.8)^(2.5-2), # t = 2.5
      lambda_3 * 0.8, # t = 3
      lambda_4 * 0.8 * (0.6/0.8)^(3.5-3), # t = 3.5
      lambda_4 * 0.6, # t = 4
      lambda_4 * 0.6 * (0.6/0.8)^(5-4), # t = 5
      lambda_4 * 0.6 * (0.6/0.8)^(6-4) # t = 6
    ), tolerance = 1e-6)

  # S(t) = c for all time points
  x = c(0.8, 0.8, 0.8)
  times = 1:3
  eval_times = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5)

  out = interp(x, times, eval_times, method = "const_surv", output = "density", add_times = FALSE)
  expect_equal(out, c(0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))

  out = interp(x, times, eval_times, method = "const_dens", output = "density", add_times = FALSE)
  expect_equal(out, c(0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))

  out = interp(x, times, eval_times, method = "const_haz", output = "density", add_times = FALSE)
  expect_equal(out, c(0.2231, 0.1995, 0.1785, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), tolerance = 1e-3)

  # edge case: S(t) = 0 for some time points
  x = c(1, 0.8, 0.8, 0, 0)
  times = 0:4
  eval_times = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5)

  out = interp(x, times, eval_times, method = "const_surv", output = "density", add_times = FALSE)
  expect_equal(out, c(0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0))

  out = interp(x, times, eval_times, method = "const_dens", output = "density", add_times = FALSE)
  expect_equal(out, c(0.2, 0.2, 0.2, 0.0, 0.0, 0.8, 0.8, 0.0, 0.0, 0.0))

  out = interp(x, times, eval_times, method = "const_haz", output = "density", add_times = FALSE)
  lambda_1 = -log(0.8/1)
  lambda_2 = -log(0.8/0.8)
  lambda_3 = -log(eps/0.8)
  lambda_4 = 0 # S(t) = 0 on the left anchor
  # all deltas are 1
  expect_equal(out,
    c(lambda_1, # t = 0
      lambda_1 * 0.8^(0.5-0), # t = 0.5
      lambda_1 * 0.8, # t = 1
      lambda_2 * 0.8 * (0.8/0.8)^(1.5-1), # t = 1.5
      lambda_2 * 0.8, # t = 2
      lambda_3 * 0.8 * (eps/0.8)^(2.5-2), # t = 2.5
      lambda_3 * 0, # t = 3
      0.0, # t = 3.5
      0.0, # t = 4
      0.0 # t = 5
    ), tolerance = 1e-6)

  # extreme cases where S at anchor is lower than internal eps: density >= 0
  x = c(0.9, 0.1, 1e-13, 1e-14, 0)
  times = 1:5
  eval_times = sort(c(0.5, times, times + 0.01))
  out = interp(x, times, eval_times, method = "const_haz", output = "density", add_times = FALSE)
  expect_all_true(out >= 0)
})
