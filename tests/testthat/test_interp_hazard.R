# we only test the vector input here since it uses the C++ matrix version internally
test_that("interp() hazard works", {
  eps = 1e-12

  # S(t) all different, 3 time points, t1 == 0
  x = c(1, 0.8, 0.5)
  times = c(0, 1, 2)
  eval_times = c(0, 0.5, 1, 1.5, 2, 3)

  out = interp(x, times, eval_times, method = "const_surv", output = "hazard", add_times = FALSE)
  expect_equal(out, c(0.0, 0.0, 0.2, 0.0, 1 - 0.5/0.8, 0.0))

  out = interp(x, times, eval_times, method = "const_dens", output = "hazard", add_times = FALSE)
  f1 = (1 - 0.8) / (1 - 0)
  f2 = (0.8 - 0.5) / (2 - 1)
  expect_equal(out,
    c(f1, # t = 0
      f1 /(1 - f1 *(0.5 - 0)), # t = 0.5
      f1 / 0.8, # t = 1
      f2 /(0.8 - f2 * (1.5 - 1)), # t = 1.5
      f2 / 0.5, # t = 2
      f2 / (0.5 - f2 * (3 - 2)) # t = 3 (S(3) > 0, so we get a non-zero hazard there)
    ), tolerance = 1e-6)

  out = interp(x, times, eval_times, method = "const_haz", output = "hazard", add_times = FALSE)
  lambda_1 = -log(0.8 / 1)
  lambda_2 = -log(0.5 / 0.8)
  expect_equal(out, c(lambda_1, lambda_1, lambda_1, lambda_2, lambda_2, lambda_2),
               tolerance = 1e-6)

  # S(t) all different, 3 time points, t1 > 0
  x = c(0.9, 0.7, 0.4)
  times = c(1, 2, 3)
  eval_times = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5)

  out = interp(x, times, eval_times, method = "const_surv", output = "hazard", add_times = FALSE)
  expect_equal(out, c(0.0, 0.0, 1 - 0.9/1.0, 0.0, 1 - 0.7/0.9, 0.0, 1 - 0.4/0.7, 0.0, 0.0),
               tolerance = 1e-6)

  out = interp(x, times, eval_times, method = "const_dens", output = "hazard", add_times = FALSE)
  f1 = (1 - 0.9) / (1 - 0)
  f2 = (0.9 - 0.7) / (2 - 1)
  f3 = (0.7 - 0.4) / (3 - 2)
  expect_equal(out,
    c(f1, # t = 0
      f1 / (1 - f1 * (0.5 - 0)), # t = 0.5
      f1 / 0.9, # t = 1
      f2 / (0.9 - f2 * (1.5 - 1)), # t = 1.5
      f2 / 0.7, # t = 2
      f3 / (0.7 - f3 * (2.5 - 2)), # t = 2.5
      f3 / 0.4, # t = 3
      f3 / (0.4 - f3 * (4 - 3)), # t = 4
      0.0 # t = 5 (S(5) = 0, so hazard is 0)
    ), tolerance = 1e-6)

  out = interp(x, times, eval_times, method = "const_haz", output = "hazard", add_times = FALSE)
  lambda_1 = -log(0.9 / 1)
  lambda_2 = -log(0.7 / 0.9)
  lambda_3 = -log(0.4 / 0.7)
  expect_equal(out, c(lambda_1, lambda_1, lambda_1, lambda_2, lambda_2, lambda_3,
    lambda_3, lambda_3, lambda_3), tolerance = 1e-6)

  # S(t) all different, 2 time points, t1 > 0
  x = c(0.8, 0.3)
  times = c(1, 2)
  eval_times = c(0, 0.5, 1, 1.5, 2, 2.5, 3)

  out = interp(x, times, eval_times, method = "const_surv", output = "hazard", add_times = FALSE)
  expect_equal(out, c(0.0, 0.0, 1 - 0.8/1, 0.0, 1 - 0.3/0.8, 0.0, 0.0))

  out = interp(x, times, eval_times, method = "const_dens", output = "hazard", add_times = FALSE)
  f1 = (1 - 0.8) / (1 - 0)
  f2 = (0.8 - 0.3) / (2 - 1)
  expect_equal(out,
    c(f1, # t = 0
      f1 / (1 - f1 * (0.5 - 0)), # t = 0.5
      f1 / 0.8, # t = 1
      f2 / (0.8 - f2 * (1.5 - 1)), # t = 1.5
      f2 / 0.3, # t = 2
      f2 / (0.3 - f2 * (2.5 - 2)), # t = 2.5 (S(2.5) > 0)
      0.0 # t = 3 (S(3) = 0, so hazard is 0)
    ), tolerance = 1e-6)

  out = interp(x, times, eval_times, method = "const_haz", output = "hazard", add_times = FALSE)
  lambda_1 = -log(0.8 / 1)
  lambda_2 = -log(0.3 / 0.8)
  expect_equal(out, c(lambda_1, lambda_1, lambda_1, lambda_2, lambda_2, lambda_2, lambda_2),
               tolerance = 1e-6)

  # 1 time point only, at t = 0
  out = interp(x = 1, times = 0, eval_times = c(0, 1), method = "const_surv",
               output = "hazard", add_times = FALSE)
  expect_equal(out, c(0, 0))

  out = interp(x = 1, times = 0, eval_times = c(0, 1), method = "const_dens",
               output = "hazard", add_times = FALSE)
  expect_equal(out, c(0, 0))

  out = interp(x = 1, times = 0, eval_times = c(0, 1), method = "const_haz",
               output = "hazard", add_times = FALSE)
  expect_equal(out, c(0, 0))

  # 1 time point > 0, with S > 0
  out = interp(x = 0.4, times = 1, eval_times = c(0, 0.5, 1, 1.5, 2),
               method = "const_surv", output = "hazard", add_times = FALSE)
  expect_equal(out, c(0, 0, 1 - 0.4/1, 0, 0))

  out = interp(x = 0.4, times = 1, eval_times = c(0, 0.5, 1, 1.5, 2),
               method = "const_dens", output = "hazard", add_times = FALSE)
  f1 = (1 - 0.4) / (1 - 0)
  expect_equal(out,
    c(f1, # t = 0
      f1 / (1 - f1 * (0.5 - 0)), # t = 0.5
      f1 / 0.4, # t = 1
      f1 / (0.4 - f1 * (1.5 - 1)), # t = 1.5, S(1.5) > 0
      0.0 # S is already 0 at t = 1, so hazard is 0 at t > 1
    ), tolerance = 1e-6)

  out = interp(x = 0.4, times = 1, eval_times = c(0, 0.5, 1, 1.5, 2, 10, 20),
               method = "const_haz", output = "hazard", add_times = FALSE)
  lambda_1 = -log(0.4)
  expect_equal(out, rep(lambda_1, 7), tolerance = 1e-6)

  # 1 time point > 0, with S = 0 (super edge case)
  out = interp(x = 0, times = 1, eval_times = c(0, 0.5, 1, 1.5),
               method = "const_surv", output = "hazard", add_times = FALSE)
  expect_equal(out, c(0, 0, 1.0, 0.0))

  out = interp(x = 0, times = 1, eval_times = c(0, 0.5, 1, 1.5),
               method = "const_dens", output = "hazard", add_times = FALSE)
  f1 = (1 - 0) / (1 - 0)
  expect_equal(out, c(f1, f1 / (1 - f1 * (0.5 - 0)), 0.0, 0.0), tolerance = 1e-6)

  out = interp(x = 0, times = 1, eval_times = c(0, 0.5, 1, 1.5),
               method = "const_haz", output = "hazard", add_times = FALSE)
  # S(t) = 0 is replaced by eps for the constant hazard interpolation
  lambda_1 = -log(eps)
  expect_equal(out, c(lambda_1, lambda_1, lambda_1, 0.0), tolerance = 1e-6)

  # S(t) with tied values
  x = c(1, 0.8, 0.8, 0.6)
  times = 1:4
  eval_times = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6)

  out = interp(x, times, eval_times, method = "const_surv", output = "hazard", add_times = FALSE)
  expect_equal(out, c(0.0, 0.0, 0.0, 0.0, 1 - 0.8/1.0, 0.0, 0.0, 0.0, 1 - 0.6/0.8, 0.0, 0.0))

  out = interp(x, times, eval_times, method = "const_dens", output = "hazard", add_times = FALSE)
  f1 = (1 - 1) / (1 - 0)
  f2 = (1 - 0.8) / (2 - 1)
  f3 = (0.8 - 0.8) / (3 - 2)
  f4 = (0.8 - 0.6) / (4 - 3)
  expect_equal(out,
    c(f1, # t = 0
      f1 / (1.0 - f1 * (0.5 - 0)), # t = 0.5
      f1 / 1.0, # t = 1
      f2 / (1.0 - f2 * (1.5 - 1)), # t = 1.5
      f2 / 0.8, # t = 2
      f3 / (0.8 - f3 * (2.5 - 2)), # t = 2.5
      f3 / 0.8, # t = 3
      f4 / (0.8 - f4 * (3.5 - 3)), # t = 3.5
      f4 / 0.6, # t = 4
      f4 / (0.6 - f4 * (5 - 4)), # t = 5, S(5) > 0
      f4 / (0.6 - f4 * (6 - 4)) # t = 6, S(6) > 0
    ), tolerance = 1e-6)

  out = interp(x, times, eval_times, method = "const_haz", output = "hazard", add_times = FALSE)
  lambda_1 = -log(1 / 1)
  lambda_2 = -log(0.8 / 1)
  lambda_3 = -log(0.8 / 0.8)
  lambda_4 = -log(0.6 / 0.8)
  expect_equal(out, c(lambda_1, lambda_1, lambda_1, lambda_2, lambda_2, lambda_3, lambda_3,
                      lambda_4, lambda_4, lambda_4, lambda_4), tolerance = 1e-6)

  # S(t) = c for all time points
  x = c(0.8, 0.8, 0.8)
  times = 1:3
  eval_times = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5)

  out = interp(x, times, eval_times, method = "const_surv", output = "hazard", add_times = FALSE)
  expect_equal(out, c(0.0, 0.0, 1 - 0.8/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))

  out = interp(x, times, eval_times, method = "const_dens", output = "hazard", add_times = FALSE)
  f1 = (1 - 0.8) / (1 - 0) # no hazard after t = 1
  expect_equal(out, c(f1, f1 / (1.0 - f1 * (0.5 - 0)), f1 / 0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
               tolerance = 1e-6)

  out = interp(x, times, eval_times, method = "const_haz", output = "hazard", add_times = FALSE)
  lambda_1 = -log(0.8 / 1)
  lambda_2 = -log(0.8 / 0.8) # lambda_3 == lambda_2 since S is constant
  expect_equal(out, c(lambda_1, lambda_1, lambda_1, lambda_2, lambda_2, lambda_2,
                      lambda_2, lambda_2, lambda_2), tolerance = 1e-6)

  # edge case: S(t) = 0 for some time points
  x = c(1, 0.8, 0.8, 0, 0)
  times = 0:4
  eval_times = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5)

  out = interp(x, times, eval_times, method = "const_surv", output = "hazard", add_times = FALSE)
  expect_equal(out,
    c(0.0, # t = 0
      0.0, # t = 0.5
      1 - 0.8/1.0, # t = 1
      0.0, # t = 1.5
      1 - 0.8/0.8, # t = 2
      0.0, # t = 2.5
      1 - 0.0/1.0, # t = 3
      0.0, # t = 3.5, S == 0 from here on, so hazard is 0
      0.0, # t = 4
      0.0 # t = 5
    ), tolerance = 1e-6)

  out = interp(x, times, eval_times, method = "const_dens", output = "hazard", add_times = FALSE)
  f1 = (1 - 0.8) / (1 - 0)
  f2 = (0.8 - 0.8) / (2 - 1)
  f3 = (0.8 - 0) / (3 - 2)
  expect_equal(out,
    c(f1, # t = 0
      f1 / (1 - f1 * (0.5 - 0)), # t = 0.5
      f1 / 0.8, # t = 1
      f2 / (0.8 - f2 * (1.5 - 1)), # t = 1.5
      f2 / 0.8, # t = 2
      f3 / (0.8 - f3 * (2.5 - 2)), # t = 2.5
      0.0, # t = 3, S(3) = 0, so hazard is 0 at t >= 3
      0.0, # t = 3.5
      0.0, # t = 4
      0.0 # t = 5
    ), tolerance = 1e-6)

  out = interp(x, times, eval_times, method = "const_haz", output = "hazard", add_times = FALSE)
  lambda_1 = -log(0.8 / 1)
  lambda_2 = -log(0.8 / 0.8)
  lambda_3 = -log(eps / 0.8)
  expect_equal(out, c(lambda_1, lambda_1, lambda_1, lambda_2, lambda_2, lambda_3, lambda_3,
                      0.0, 0.0, 0.0), tolerance = 1e-6)

  # extreme cases where S at anchor is lower than internal eps: hazard >= 0
  x = c(0.9, 0.1, 1e-13, 1e-14, 0)
  times = 1:5
  eval_times = sort(c(0.5, times, times + 0.01))
  out = interp(x, times, eval_times, method = "const_haz", output = "hazard", add_times = FALSE)
  expect_all_true(out >= 0)
})
