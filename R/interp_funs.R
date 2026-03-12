#' @title Interpolate Survival Curves
#'
#' @description
#' Interpolates survival curves (vector or matrix) at new time points using
#' internal C++ interpolation functions.
#' Input must always be \emph{survival probabilities}.
#' Output can be the survival, cumulative distribution, or density functions, as well as
#' the hazard or cumulative hazard functions.
#'
#' @param x (`numeric()` | `matrix()`)\cr
#'   Survival vector or matrix (rows = observations, columns = time points).
#' @param times (`numeric()` | `NULL`)\cr
#'   Anchor time points. If `NULL`, extracted from names/colnames.
#' @param output (`character(1)`)\cr
#'   Output type: `"surv"`, `"cdf"`, `"cumhaz"`, `"density"` or `"hazard"`.
#' @param add_times (`logical(1)`)\cr
#'   If `TRUE`, attach `eval_times` as names/colnames.
#'
#' @template param_method
#' @template param_eval_times
#' @template param_check
#' @template param_eps
#' @template param_trim_duplicates
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
#' # constant S(t) interpolation
#' interp(x, times, eval_times)
#' # linear S(t) interpolation
#' interp(x, times, eval_times, method = "linear_surv")
#' # exponential S(t) interpolation (same as `method = "const_haz"`)
#' interp(x, times, eval_times, method = "exp_surv")
#'
#' # cumulative distribution with linear S(t) interpolation
#' interp(x, times, eval_times, method = "linear_surv", output = "cdf")
#'
#' # H(t) with linear interpolation
#' interp(x, times, eval_times, method = "linear_surv", output = "cumhaz")
#'
#' @export
interp = function(x,
                  times = NULL,
                  eval_times = NULL,
                  method = "const_surv",
                  output = "surv",
                  add_times = TRUE,
                  check = TRUE,
                  eps = 1e-12,
                  trim_duplicates = FALSE) {
  # quick assertions
  method = map_interp_method(method) # const_* aliases
  output = assert_choice(output, c("surv", "cdf", "cumhaz", "density", "hazard"))
  assert_flag(add_times)
  assert_flag(check)
  assert_flag(trim_duplicates)
  eval_times = assert_numeric(
    eval_times, lower = 0, unique = TRUE, sorted = TRUE,
    null.ok = TRUE, any.missing = FALSE, min.len = 1
  )
  is_mat = is.matrix(x)

  # remove flat S(t) segments
  if (trim_duplicates) {
    trimmed = trim_duplicates(x, times = times)
    x = trimmed$x
    times = trimmed$times
  }

  # optional S(t) check
  if (check) {
    times = assert_prob(x, times, type = "surv")
  } else {
    times = extract_times(x, times)
  }

  # Case: no interpolation requested => use anchor times
  if (is.null(eval_times)) {
    eval_times = times
  }
  x_mat = if (is_mat) x else matrix(x, nrow = 1)

  if (output %in% c("surv", "cdf", "cumhaz")) {
    res = if (identical(eval_times, times)) {
      x_mat # we have S(t) at the anchors already
    } else {
      c_interp_surv_mat(x_mat, times, eval_times, method)
    }
  } else if (output == "density") {
    res = c_interp_density_mat(x_mat, times, eval_times, method)
  } else if (output == "hazard") {
    res = c_interp_hazard_mat(x_mat, times, eval_times, method)
  }

  # if input was a vector, return a vector
  if (!is_mat) res = res[1, ]

  # transform S(t) => F(t) or H(t) if needed and attach time labels
  process_output(res, eval_times, output, add_times, eps)
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
