#' Interpolate Survival Curves
#'
#' Interpolates survival curves (vector or matrix) at new time points using
#' internal C interpolation functions.
#' Output can be the survival, cumulative distribution, or density functions, as well as
#' the hazard or cumulative hazard function.
#' Input must always represent survival probabilities.
#'
#' @param x (`numeric()` | `matrix()`)\cr
#'   Survival vector or matrix (rows = observations, columns = time points).
#' @param times (`numeric()` | `NULL`)\cr
#'   Original time points. If `NULL`, extracted from names/colnames.
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
#'
#' @return A numeric vector or matrix of interpolated values.
#' 
#' @examples
#' x = matrix(c(1, 0.8, 0.6,
#'              1, 0.7, 0.4),
#'            nrow = 2, byrow = TRUE)
#' times = c(0, 10, 20)
#' eval_times = c(5, 15, 25)
#' interp(x, times, eval_times, output = "surv")
#' interp(x, times, eval_times, method = "linear_surv", output = "surv")
#' @export
interp = function(x,
                  times = NULL,
                  eval_times,
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
  # access names/colnames in a generic way
  is_mat = is.matrix(x)
  name_attr = if (is_mat) colnames else names

  # optional S(t) check
  if (check) {
    times = assert_prob(x, times, type = "surv")
  } else {
    times = extract_times(x, times)
  }

  # Case: no interpolation requested
  if (is.null(eval_times)) {
    if (add_times) {
      if (is.null(name_attr(x))) {
        name_attr(x) = as.character(times)
      }
    }
    return(x)
  }

  # call C++ interpolation
  if (is_mat) {
    res = c_interp_surv_mat(x, times, eval_times, method)
  } else {
    res = c_interp_surv_vec(x, times, eval_times, method)
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
    name_attr(res) = as.character(eval_times)
  }

  res
}
