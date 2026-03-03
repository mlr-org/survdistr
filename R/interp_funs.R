#' @title Interpolate Survival Curves
#'
#' @description
#' Interpolates survival curves (vector or matrix) at new time points using
#' internal C interpolation functions.
#' Output can be the survival, cumulative distribution, or density functions, as well as
#' the hazard or cumulative hazard functions.
#' Input must always be survival probabilities.
#'
#' @param x (`numeric()` | `matrix()`)\cr
#'   Survival vector or matrix (rows = observations, columns = time points).
#' @param times (`numeric()` | `NULL`)\cr
#'   Anchor time points. If `NULL`, extracted from names/colnames.
#' @param method (`character(1)`)\cr
#'   Interpolation method: `"const_surv"` or `"linear_surv"`.
#' @param output (`character(1)`)\cr
#'   Output type: `"surv"`, `"cdf"`, or `"cumhaz"`.
#' @param add_times (`logical(1)`)\cr
#'   If `TRUE`, attach `eval_times` as names/colnames.
#'
#' @template param_eval_times
#' @template param_check
#' @template param_eps
#' @templateVar eps 1e-6
#'
#' @return A numeric vector or matrix of interpolated values.
#' 
#' @examples
#' x = matrix(c(1, 0.8, 0.6,
#'              1, 0.7, 0.4),
#'            nrow = 2, byrow = TRUE)
#' times = c(0, 10, 20)
#' eval_times = c(5, 15, 25)
#'
#' # S(t) with constant interpolation
#' interp(x, times, eval_times)
#' # S(t) with linear interpolation
#' interp(x, times, eval_times, method = "linear_surv")
#' # H(t) with linear interpolation
#' interp(x, times, eval_times, method = "linear_surv", output = "cumhaz")
#' @export
interp = function(x,
                  times = NULL,
                  eval_times = NULL,
                  method = "const_surv",
                  output = "surv",
                  add_times = TRUE,
                  check = TRUE,
                  eps = 1e-6) {
  # quick assertions
  method = assert_choice(method, c("const_surv", "linear_surv"))
  output = assert_choice(output, c("surv", "cdf", "cumhaz"))
  assert_flag(add_times)
  assert_flag(check)
  eval_times = assert_numeric(
    eval_times, lower = 0, unique = TRUE, sorted = TRUE,
    null.ok = TRUE, any.missing = FALSE, min.len = 1
  )
  is_mat = is.matrix(x)

  # optional S(t) check
  if (check) {
    times = assert_prob(x, times, type = "surv")
  } else {
    times = extract_times(x, times)
  }

  # Case: no interpolation requested
  if (is.null(eval_times)) {
    if (add_times) {
      if (is_mat && is.null(colnames(x))) colnames(x) = as.character(times)
      if (!is_mat && is.null(names(x))) names(x) = as.character(times)
    }
    return(x)
  }

  # call C++ interpolation
  if (is_mat) {
    res = c_interp_surv_mat(x, times, eval_times, method)
  } else {
    res = c_interp_surv_mat(matrix(x, nrow = 1), times, eval_times, method)[1, ]
  }

  # transform output
  if (output == "cdf") {
    res = 1 - res
  } else if (output == "cumhaz") {
    # avoid -Inf for zero survival probabilities
    res = -log(pmax(res, eps))
  }

   # attach time labels
  if (add_times) {
    if (is_mat) {
      colnames(res) = as.character(eval_times)
    } else {
      names(res) = as.character(eval_times)
    }
  }

  res
}

#' Interpolate CIF matrix
#'
#' Interpolates cumulative incidence (CIF) functions (corresponding to one competing event only)
#' using left-continuous constant interpolation.
#'
#' @param x (`matrix()`)\cr
#'   CIF matrix (rows = observations, columns = time points).
#' @param times (`numeric()` | `NULL`)\cr
#'   Anchor time points. If `NULL`, extracted from `colnames(x)`.
#' @param add_times (`logical(1)`)\cr
#'   If `TRUE`, attach `eval_times` as colnames in the output matrix.
#' @template param_eval_times
#' @template param_check
#'
#' @return Interpolated CIF matrix.
#' @export
interp_cif = function(x, times = NULL, eval_times = NULL, add_times = TRUE, check = TRUE) {
  # quick assertions
  assert_flag(add_times)
  assert_flag(check)
  eval_times = assert_numeric(
    eval_times, lower = 0, unique = TRUE, sorted = TRUE,
    null.ok = TRUE, any.missing = FALSE, min.len = 1
  )

  # optional CIF(t) check
  if (check) {
    times = assert_prob_matrix(x, times, type = "cif")
  } else {
    times = extract_times(x, times)
  }

  # Case: no interpolation requested
  if (is.null(eval_times)) {
    if (add_times) {
      if (is.null(colnames(x))) {
        colnames(x) = as.character(times)
      }
    }
    return(x)
  }

  # call C++ interpolation
  res = c_interp_cif_mat(x, times, eval_times)

  if (add_times) {
    colnames(res) = as.character(eval_times)
  }

  res
}
