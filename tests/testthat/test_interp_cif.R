test_that("interp_cif works", {
  # Simple CIF matrix: 2 observations, 3 time points
  cif = matrix(c(0, 0.2, 0.5,
                 0, 0.3, 0.6),
                nrow = 2, byrow = TRUE)
  times = c(1, 5, 10)
  eval_times = c(0, 2, 5, 8, 12)
  # t < 1 → 0
  # 1 ≤ t < 5 → CIF(1) = 0 (row1), 0 (row2)
  # 5 ≤ t < 10 → CIF(5) = 0.2 (row1), 0.3 (row2)
  # t ≥ 10 → CIF(10) = 0.5 (row1), 0.6 (row2)

  expected = matrix(c(0, 0, 0.2, 0.2, 0.5,
                      0, 0, 0.3, 0.3, 0.6),
                     nrow = 2, byrow = TRUE)
  colnames(expected) = as.character(eval_times)

  res = interp_cif(cif, times, eval_times, add_times = TRUE, check = FALSE)
  expect_equal(res, expected)
})

test_that("interp_cif handles boundary conditions correctly", {
  cif = matrix(c(0, 0.4, 0.7), nrow = 1)  # one observation
  times = c(2, 5, 8)
  eval_times = c(0, 2, 3, 5, 7, 8, 10)

  expected = matrix(c(0, 0, 0, 0.4, 0.4, 0.7, 0.7), nrow = 1)
  colnames(expected) = as.character(eval_times)

  res = interp_cif(cif, times, eval_times, add_times = TRUE, check = FALSE)
  expect_equal(res, expected)
})

test_that("interp_cif handles single time point", {
  cif = matrix(c(0, 0.2, 0.5), nrow = 3) # 3 observations
  times = 5 # 1 time point
  eval_times = c(0, 2, 5, 7)

  # Expected: t < 5 → 0, t ≥ 5 → value
  expected = matrix(c(0, 0, 0, 0,
                      0, 0, 0.2, 0.2,
                      0, 0, 0.5, 0.5), nrow = 3, byrow = TRUE)

  res = interp_cif(cif, times, eval_times, check = FALSE)
  expect_equal(res, expected, ignore_attr = TRUE)
})


test_that("interp_cif returns original matrix when eval_times is NULL", {
  cif = gen_random_mat(3, 4, type = "cif")
  times = c(1, 3, 5, 7)

  # Without adding times
  res = interp_cif(cif, times, eval_times = NULL, add_times = FALSE)
  expect_equal(res, cif)
  expect_null(colnames(res))

  # With adding times (if colnames missing)
  res2 = interp_cif(cif, times, eval_times = NULL, add_times = TRUE)
  expect_equal(res2, cif, ignore_attr = TRUE)  # values unchanged
  expect_equal(colnames(res2), as.character(times))
})

test_that("interp_cif extracts times from colnames if not provided", {
  cif = gen_random_mat(2, 3, type = "cif")
  times_char = c("1", "5", "10")
  colnames(cif) = times_char

  eval_times = c(0, 2, 5, 8, 12)
  # Use check=FALSE to skip validation (so times are extracted from colnames)
  res = interp_cif(cif, times = NULL, eval_times, check = FALSE)

  # Compute expected manually or compare with providing times explicitly
  res_explicit = interp_cif(cif, times = as.numeric(times_char), eval_times, check = FALSE)
  expect_equal(res, res_explicit)
})

test_that("interp_cif checks", {
  # Non-decreasing violation
  bad_cif = matrix(c(0, 0.5, 0.3), nrow = 1)
  times = c(1, 2, 3)
  expect_error(interp_cif(bad_cif, times, eval_times = c(0, 1.5), check = TRUE))

  # Values outside [0,1]
  bad_cif2 = matrix(c(-0.1, 0.2, 0.8), nrow = 1)
  expect_error(interp_cif(bad_cif2, times, eval_times = c(0, 1.5), check = TRUE))
})
