test_that("convert_to_surv() works", {
  times = c(1, 3, 6)

  # discrete density
  disc_dens = c(0.1, 0.2, 0.15)
  expect_equal(
    convert_to_surv(disc_dens, times = times, input = "disc_dens"),
    c(0.9, 0.7, 0.55)
  )

  # continuous density (Delta = 1, 2, 3)
  cont_dens = c(0.1, 0.05, 0.1)
  expect_equal(
    convert_to_surv(cont_dens, times = times, input = "cont_dens"),
    c(0.9, 0.8, 0.5)
  )

  # discrete hazard
  disc_haz = c(0.1, 0.2, 0.5)
  expect_equal(
    convert_to_surv(disc_haz, times = times, input = "disc_haz"),
    c(0.9, 0.72, 0.36)
  )

  # continuous hazard (delta = 1, 2, 3)
  cont_haz = c(0.1, 0.05, 0.2)
  expect_equal(
    convert_to_surv(cont_haz, times = times, input = "cont_haz"),
    exp(-c(0.1, 0.2, 0.8))
  )

  # matrix input works
  mat = matrix(c(0.1, 0.2, 0.2, 0.3), nrow = 2, byrow = TRUE)
  res = convert_to_surv(mat, times = c(1, 2), input = "disc_haz")
  expect_equal(dim(res), dim(mat))
  expect_equal(res, matrix(c(0.9, 0.72, 0.8, 0.56), nrow = 2, byrow = TRUE))

  # non-proper times gives an error always
  expect_error(convert_to_surv(c(0.1, 0.12), times = c(1, -2)), "is not >= 0")

  # several checks and corrections
  expect_silent(convert_to_surv(c(0.1, -0.1), times = c(1, 2), input = "disc_dens", check = FALSE))
  expect_error(convert_to_surv(c(0.1, -0.1), times = c(1, 2), input = "disc_dens", check = TRUE),
               "Density probabilities must be in [0,1].", fixed = TRUE)
  res = expect_warning(convert_to_surv(c(0.1, -0.2), times = c(1, 2), input = "disc_dens", check = FALSE,
                       clamp_surv = TRUE), "Survival values outside")
  expect_silent(convert_to_surv(c(4, 3), times = c(1, 2), input = "disc_haz", check = FALSE))
  expect_error(convert_to_surv(c(4, 3), times = c(1, 2), input = "disc_haz", check = TRUE),
               "Hazard probabilities must be in [0,1].", fixed = TRUE)
  expect_equal(convert_to_surv(c(0.8, 0.6, 0.8, 0.8), times = 1:4, input = 'disc_haz',
               clamp_surv = TRUE, eps = 0.01), c(0.2, 0.08, 0.016, 0.01))

  # S(0) = 1 for discrete inputs errors
  expect_silent(convert_to_surv(c(0.1, 0.2), times = c(0, 1), input = "disc_dens", check = FALSE))
  expect_error(convert_to_surv(c(0.1, 0.2), times = c(0, 1), input = "disc_dens", check = TRUE),
               "times must be strictly positive")
  expect_silent(convert_to_surv(c(0.1, 0.2), times = c(0, 1), input = "disc_haz", check = FALSE))
  expect_error(convert_to_surv(c(0.1, 0.2), times = c(0, 1), input = "disc_haz", check = TRUE),
               "times must be strictly positive")
  # S(0) = 1 for continuous inputs works
  expect_equal(convert_to_surv(c(0.1, 0.2), times = c(0, 1), input = "cont_dens")[1], 1)
  expect_equal(convert_to_surv(c(0.1, 0.2), times = c(0, 1), input = "cont_haz")[1], 1)
})
