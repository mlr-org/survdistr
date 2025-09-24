#' Interpolate Survival / CDF / CIF Matrices
#'
#' Wrapper around the internal C++ interpolation function \code{c_mat_interp}.
#' Performs input validation before calling the underlying C++ code.
#' Can be used for survival, cumulative distribution (CDF), or cumulative
#' incidence (CIF) matrices.
#'
#' @param x (`matrix()`)\cr Survival/CDF/CIF matrix with rows as observations
#'  and columns as time points.
#' @param times (`numeric()` | `NULL`)\cr
#'  Original time points corresponding to columns of `x`.
#' @template param_eval_times
#' @template param_constant
#' @template param_type
#' @template param_add_times
#' @template param_check
#'
#' @return A numeric matrix with the same number of rows as `x` and number of
#'   columns equal to `length(eval_times)`.
#'
#' @examples
#' x = matrix(c(1, 0.8, 0.6,
#'              1, 0.7, 0.4),
#'            nrow = 2, byrow = TRUE)
#' times = c(0, 10, 20)
#' eval_times = c(5, 15, 25, 15) # duplicates & unordered
#' mat_interp(x, times, eval_times, constant = TRUE, type = "surv")
#' @export
mat_interp = function(x, times = NULL, eval_times = NULL, constant = TRUE, type = "surv",
                      add_times = TRUE, check = TRUE) {
  # quick assertions
  assert_flag(constant)
  type = assert_choice(type, c("surv", "cdf", "cif"))
  assert_flag(add_times)
  assert_flag(check)
  eval_times = assert_numeric(eval_times, lower = 0, any.missing = FALSE,
                              null.ok = TRUE, min.len = 1)

  # Optional matrix check
  if (check) {
    times = assert_prob_matrix(x, times, type = type)
  } else {
    # at least derive and check the time points
    times = extract_times(x, times)
  }

  # case: no interpolation requested
  if (is.null(eval_times)) {
    if (add_times && is.null(colnames(x))) {
      colnames(x) = as.character(times)
    }
    return(x)
  }

  # unique + sorted eval_times for C++
  eval_times_unique = sort(unique(eval_times))

  # call C++ interpolation
  mat = c_mat_interp(x, times, eval_times_unique, constant, type)

  # map back to requested order (with duplicates) if necessary
  if (!identical(eval_times, eval_times_unique)) {
    idx = match(eval_times, eval_times_unique)
    mat = mat[, idx, drop = FALSE]
  }

  if (add_times) {
    colnames(mat) = as.character(eval_times)
  }

  mat
}

#' Interpolate a Survival / CDF / CIF Vector
#'
#' Wrapper around the internal C++ interpolation function \code{c_vec_interp}.
#' Performs input validation before calling the underlying C++ code.
#' Can be used for survival, cumulative distribution (CDF), or cumulative
#' incidence (CIF) curves (vectors).
#'
#' @param x (`numeric()`)\cr
#'   Survival/CDF/CIF vector at given time points.
#'   Optionally named with the corresponding times.
#' @param times (`numeric()` | `NULL`)\cr
#'   Original time points corresponding to `x`.
#'   If `NULL`, extracted from `names(x)`.
#' @template param_eval_times
#' @template param_constant
#' @template param_type
#' @template param_add_times
#' @param check (`logical(1)`)\cr
#'   If `TRUE` (default), perform simple validation (range, monotonicity, and bounds).
#'   Set to `FALSE` to skip checks (NOT recommended for external use).
#'
#' @return A numeric vector of length `length(eval_times)` with interpolated values.
#'
#' @examples
#' surv_vec = c(1, 0.8, 0.6)
#' names(surv_vec) = c(0, 10, 20)
#' eval_times = c(5, 15, 25)
#' vec_interp(surv_vec, eval_times = eval_times, type = "surv")
#' @export
vec_interp = function(x, times = NULL, eval_times = NULL, constant = TRUE,
                      type = "surv", add_times = TRUE, check = TRUE) {
  # Quick assertions
  assert_numeric(x, any.missing = FALSE, min.len = 1)
  assert_flag(constant)
  type = assert_choice(type, c("surv", "cdf", "cif"))
  assert_flag(add_times)
  assert_flag(check)
  eval_times = assert_numeric(eval_times, lower = 0, any.missing = FALSE,
                              null.ok = TRUE, min.len = 1)

  # get `times` from argument or names(x)
  if (is.null(times)) {
    if (is.null(names(x))) {
      stop("If 'times' is NULL, names(x) must provide the time points.")
    }
    times = assert_numeric(as.numeric(names(x)),
                           lower = 0, unique = TRUE, sorted = TRUE,
                           any.missing = FALSE)
  } else {
    times = assert_numeric(times,
                           lower = 0, unique = TRUE, sorted = TRUE,
                           any.missing = FALSE, null.ok = FALSE)
    if (length(x) != length(times)) {
      stop("Length of 'x' must match length of 'times'.")
    }
  }

  # Simple check (if requested)
  if (check) {
    if (any(x < 0 | x > 1)) stop("Values must be in [0,1].")
    diffs = diff(x)
    if (type == "surv" && any(diffs > 0)) stop("Survival must be non-increasing.")
    if (type %in% c("cdf", "cif") && any(diffs < 0)) stop("CDF/CIF must be non-decreasing.")
    if (times[1] == 0) {
      if (type == "surv" && x[1] != 1) stop("S(0) must equal 1.")
      if (type %in% c("cdf", "cif") && x[1] != 0) stop("CDF/CIF(0) must equal 0.")
    }
  }

  # Case: no interpolation requested
  if (is.null(eval_times)) {
    if (add_times && is.null(names(x))) {
      names(x) = as.character(times)
    }
    return(x)
  }

  # Unique + sorted eval_times for C++
  eval_times_unique = sort(unique(eval_times))

  # call C++ interpolation
  vec = c_vec_interp(x, times, eval_times_unique, constant, type)

  # Map back to requested eval_times
  if (!identical(eval_times, eval_times_unique)) {
    idx = match(eval_times, eval_times_unique)
    vec = vec[idx]
  }

  if (add_times) {
    names(vec) = as.character(eval_times)
  }

  vec
}
